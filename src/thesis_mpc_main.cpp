#include <stddef.h>
#include <stdio.h>    
#include <math.h> 
#include <cmath>
#include <sstream>             
#include <iostream>
#include <fstream>
#include <tf2/LinearMath/Quaternion.h>
// #include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseArray.h>
#include <vector>
//mpc controller
#include "../include/shao_thesis_mpc/thesis_mpc.h"               
//ROS
#include "ros/ros.h"
#include "std_msgs/String.h"
//Surscribe
#include <nav_msgs/Odometry.h>
#include "sensor_msgs/Imu.h"
#include <mavros_msgs/State.h>
//Publish
#include <mavros_msgs/AttitudeTarget.h>
#include <mavros_msgs/ParamSet.h>
#include <mavros_msgs/ParamValue.h>

// PPO controller
// #include <onnxruntime/core/providers/cpu/cpu_provider_factory.h>
// #include <onnxruntime_cxx_api.h>

using namespace std;

constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }
double pix_Y2mpc_Y(double x) {
    if (x - pi() / 2 < -pi()) {
        return -(x - pi() / 2 + 2 * pi());
    } else {
        return -(x - pi() / 2);
    }
}
double norm(double x, double Finite) {
    if ((x - (-Finite)) / (Finite - (-Finite)) * 2 - 1 > 1) {
        return 1;
    } else if ((x - (-Finite)) / (Finite - (-Finite)) * 2 - 1 < -1) {
        return -1;
    } else {
        return (x - (-Finite)) / (Finite - (-Finite)) * 2 - 1;
    }
}


double calculateRate(double currentValue, double previousValue) {
    return (currentValue - previousValue) / 0.1;
}

struct Point {
    double x, y;
};

class Listener {
public:
    double qua_w, qua_x, qua_y, qua_z;
    double siny_cosp, cosy_cosp, sinp;
    double UAV_yaw, UAV_phi;
    double UAV_E, UAV_N, UAV_U, UAV_p, UAV_q, UAV_r;
    string guided_state;
    vector<Point> path;

    void imu_callback(const sensor_msgs::Imu::ConstPtr& msg);
    void pos_callback(const nav_msgs::Odometry::ConstPtr& msg);
    void guided_callback(const mavros_msgs::State::ConstPtr& msg);
    void trajectory_callback(const geometry_msgs::PoseArray::ConstPtr& msg);
};

void Listener::imu_callback(const sensor_msgs::Imu::ConstPtr& msg) {
    std_msgs::String pub_str;
    std::stringstream ss;
    ss << "controller heard: x: " << msg->orientation.x << " y: " << msg->orientation.y << " z: " << msg->orientation.z << " w: " << msg->orientation.w;
    pub_str.data = ss.str();

    qua_w = msg->orientation.w;
    qua_x = msg->orientation.x;
    qua_y = msg->orientation.y;
    qua_z = msg->orientation.z; 

    UAV_p = msg->angular_velocity.x;
    UAV_q = msg->angular_velocity.y;
    UAV_r = msg->angular_velocity.z;

    siny_cosp = 2 * (qua_w * qua_z + qua_x * qua_y);
    cosy_cosp = 1 - 2 * (qua_y * qua_y + qua_z * qua_z);
    UAV_yaw = std::atan2(siny_cosp, cosy_cosp);
    sinp = 2 * (qua_w * qua_y - qua_z * qua_x);
    if (std::abs(sinp) >= 1) // 防止因數值超出範圍導致的誤差
        UAV_phi = std::copysign(M_PI / 2, sinp); // 使用90度或-90度
    else
        UAV_phi = std::asin(sinp);
}

void Listener::pos_callback(const nav_msgs::Odometry::ConstPtr& msg) {
    UAV_E = msg->pose.pose.position.x;   // East(x)
    UAV_N = msg->pose.pose.position.y;   // North(y)
    UAV_U = msg->pose.pose.position.z;   // Up(z)
}

void Listener::guided_callback(const mavros_msgs::State::ConstPtr& msg) {
    guided_state = msg->mode;
}

void Listener::trajectory_callback(const geometry_msgs::PoseArray::ConstPtr& msg) {
    path.clear(); // Clear previous path points

    for (const auto& pose : msg->poses) {
        Point new_point;
        new_point.x = pose.position.x;
        new_point.y = pose.position.y;
        path.push_back(new_point);
    }
}


std::vector<double> crossProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
    return {
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0]
    };
}

std::tuple<double, double, double, double, double, double, double, int> fcn(
    double P_N, double P_E, const std::vector<Point>& path,
    double psi, double uav_r)
{
    double v = 1.5; // Set speed

    double uav_k = uav_r / 15;
    int numWaypoints = path.size(); // Number of points in the path

    // Calculate the distances from the current position to all waypoints
    std::vector<double> distances(numWaypoints);
    for (int i = 0; i < numWaypoints; ++i) {
        distances[i] = std::sqrt(std::pow(path[i].y - P_N, 2) + std::pow(path[i].x - P_E, 2));
    }

    // Find the index of the closest waypoint
    int current_index = std::min_element(distances.begin(), distances.end()) - distances.begin();

    // Calculate the target point index
    std::vector<double> distances_after_nearest(distances.begin() + current_index + 1, distances.end());
    double target_distance = 80;
    std::vector<double> diff_to_target_after_nearest(distances_after_nearest.size());
    std::transform(distances_after_nearest.begin(), distances_after_nearest.end(), diff_to_target_after_nearest.begin(),
                   [target_distance](double d) { return std::abs(d - target_distance); });

    int index_after_nearest = std::min_element(diff_to_target_after_nearest.begin(), diff_to_target_after_nearest.end()) - diff_to_target_after_nearest.begin();
    int i = index_after_nearest + current_index + 1;

    double N = path[i].y; // Predicted north position
    double E = path[i].x; // Predicted east position

    // Calculate cross-track error
    double x1 = path[current_index].x;
    double y1 = path[current_index].y;
    double x2 = path[current_index + 1].x;
    double y2 = path[current_index + 1].y;

    double D = std::sqrt(std::pow(P_N - y1, 2) + std::pow(P_E - x1, 2));
    double psi_r = std::atan2(y2 - y1, x2 - x1);
    double psi_c = std::atan2(P_N - y1, P_E - x1);
  

    double d = std::abs(D * std::sin(psi_r - psi_c));

    // Calculate curvature using three points
    double x1_C = path[current_index + 1].x;
    double y1_C = path[current_index + 1].y;
    double x2_C = path[current_index + 2].x;
    double y2_C = path[current_index + 2].y;
    double x3_C = path[current_index + 3].x;
    double y3_C = path[current_index + 3].y;

    double a = std::sqrt(std::pow(x2_C - x3_C, 2) + std::pow(y2_C - y3_C, 2));
    double b = std::sqrt(std::pow(x1_C - x3_C, 2) + std::pow(y1_C - y3_C, 2));
    double c = std::sqrt(std::pow(x1_C - x2_C, 2) + std::pow(y1_C - y2_C, 2));

    double curvature;
    double r;
    if ((x1_C == x2_C && x2_C == x3_C) || (y1_C == y2_C && y2_C == y3_C)) {
        curvature = 0;
        r = 0;
    } else {
        double cosB = (b * b + c * c - a * a) / (2 * b * c);
        std::vector<double> AB = {x2_C - x1_C, y2_C - y1_C, 0};
        std::vector<double> AC = {x3_C - x1_C, y3_C - y1_C, 0};
        std::vector<double> cross = crossProduct(AB, AC);

        if (std::abs(cosB) > 1) {
            curvature = 0;
            r = 0;
        } else {
            double sinB = std::sqrt(1 - cosB * cosB);
            if (sinB > 1e-6) {
                r = a / (2 * sinB);
                curvature = 1 / r;
                curvature *= (cross[2] >= 0) ? 1 : -1;
            } else {
                curvature = 0;
                r = 0;
            }
        }
    }

    r = std::abs(r);

    // Calculate heading angle change and other control variables
    double psi_cmd = std::atan2(N - P_N, E - P_E);
    double L1 = std::sqrt(std::pow(N - P_N, 2) + std::pow(E - P_E, 2));
    double a_cmd = 2 * std::pow(v, 2) / L1 * std::sin(psi_cmd - psi);
    double phi = std::atan2(-a_cmd, 9.81) * (180.0 / M_PI);

    double heta = psi_cmd - psi;
    double R = L1 / (2 * std::sin(heta));
    double delta_psi = (a_cmd / v) / 0.1;

    return {N, E, d, delta_psi, heta, uav_k, curvature, current_index};
}

thesis_mpc rtObj;

// Guidanse Law Parameters
double path_psi = 0;
double path_N_start = 0;
double path_E_start = 0;
double over_circle_N = 0;
double over_circle_E = 0;
double path_N_mid;
double path_E_mid;
double rotation_heading;
double radius = 100;
double gostraight = 30;
// double L1 = 70;
// double V = 15;
double Targetpoint_1[2]; // [North,East]
double Targetpoint_2[2];
double path_N_end;
double path_E_end;
double line_a;
double line_b;
double line_c;
double Inter_d;
double Inter_A;
double Inter_h;
double Inter_N;
double Inter_E;
double choose_Tp1;
double choose_Tp2;
double ahead_uav_E;
double ahead_uav_N;
double a_cmd;
double distance_error;
double r_cmd;
double UAV_roll_cmd;
double UAV_yaw_cmd;
double distance_error_past = 0;

// Atitude Hold Parameter
double Up_cmd;
double Lon_P;
double UAV_pitch_cmd;

double r_e;
float N_r_e;
float N_r;
float N_d;
float N_d_dot;
float N_N;
float N_E;
float N_L1;
float N_heta;
float N_R;
float N_phi;
float N_phi_r;
float phi_r;

double Mean;
double standarddeviation;
bool modechange = false;

int main(int argc, char** argv) {

    ROS_INFO("A");

    ros::init(argc, argv, "thesis_main");
    ros::NodeHandle n;
    Listener listener;
    ros::Subscriber sub_imu = n.subscribe<sensor_msgs::Imu>("mavros/imu/data", 1, &Listener::imu_callback, &listener);
    ros::Subscriber sub_position = n.subscribe<nav_msgs::Odometry>("mavros/global_position/local", 1, &Listener::pos_callback, &listener);
    ros::Subscriber sub_guided = n.subscribe<mavros_msgs::State>("mavros/state", 1, &Listener::guided_callback, &listener);
    ros::Subscriber sub_trajectory = n.subscribe("/predicted_trajectory", 1, &Listener::trajectory_callback, &listener);
    
    //what is this??
    mavros_msgs::AttitudeTarget r;
    mavros_msgs::ParamSet s;
    mavros_msgs::ParamSet e;

    ros::Publisher pub = n.advertise<mavros_msgs::AttitudeTarget>("mavros/setpoint_raw/attitude", 1);
    ros::ServiceClient pubth = n.serviceClient<mavros_msgs::ParamSet>("mavros/param/set");
    ros::Rate rate(10.0);

    while (ros::ok()) {

        if (listener.guided_state == "GUIDED") {
            s.request.param_id = "SERVO3_MIN";
            s.request.value.integer = 1600;

            e.request.param_id = "SERVO4_FUNCTION";
            e.request.value.integer = 0;
            r.header.stamp = ros::Time::now();

            double P_E = listener.UAV_E, P_N = listener.UAV_N; // 當前位置

            double psi =listener.UAV_yaw; // 當前航向角
            double uav_r = listener.UAV_r;
            
            auto result = fcn(listener.UAV_N, listener.UAV_E, listener.path, listener.UAV_yaw, listener.UAV_r);

            // 提取元組中的值
            double N = std::get<0>(result);
            double E = std::get<1>(result);
            double d = std::get<2>(result);
            double delta_psi = std::get<3>(result);
            double heta = std::get<4>(result);
            double uav_k = std::get<5>(result);
            double curvature = std::get<6>(result);
            int current_index = std::get<7>(result);

            double L1 = sqrt(pow(N - P_N, 2) + pow(E - P_E, 2));
            double a_cmd = 2 * pow(1.5, 2) / L1 * sin(heta);
            // double UAV_roll_cmd = atan2(-a_cmd, 9.81) * (180.0 / M_PI);

            double uav_yaw_function;
            uav_yaw_function = pix_Y2mpc_Y(listener.UAV_yaw);
            double distance_error = d;
            double distance_error_dot = calculateRate(distance_error, distance_error_past);
            distance_error_past = distance_error;

            // 輸出結果
            std::cout << "=====================================================================" << std::endl;
            std::cout << " , a_cmd: " << a_cmd << std::endl;
            std::cout << " , delta_psi: " << delta_psi << " , Cross-track distance: " << distance_error << std::endl;
            std::cout << " ,N: " << N << " ,E: " << E << " ,heta: " << heta << std::endl;
            std::cout << " uav_N: " << P_N << " uav_E: " << P_E << " UAV_psi: " <<  rad2deg(psi)<< "UAV_yaw_function"<< rad2deg(uav_yaw_function) << std::endl;
            std::cout << "=====================================================================" << std::endl;

            // Path Following Position Controller
            rtObj.rtU.uav_r = listener.UAV_r;
            rtObj.rtU.r_cmd = delta_psi;
            rtObj.step();
            UAV_roll_cmd = rtObj.rtY.Out1;
            std::cout << "Roll_cmd: " <<rad2deg(-rtObj.rtY.Out1) <<" Roll_cmd_rad: "<<-rtObj.rtY.Out1<<  std::endl;
            // ofs  << "Roll_cmd: " <<rad2deg(-rtObj.rtY.Out1) <<" Roll_cmd_rad: "<<-rtObj.rtY.Out1<< "\n";
            UAV_yaw_cmd = 0;

            // Height Control
            Up_cmd = 70;
            Lon_P = 0.6;
            if ((Up_cmd - listener.UAV_U) * Lon_P > 10) {
                UAV_pitch_cmd = 10 * pi() / 180;
            } else if ((Up_cmd - listener.UAV_U) * Lon_P < -10) {
                UAV_pitch_cmd = -10 * pi() / 180;
            } else {
                UAV_pitch_cmd = (Up_cmd - listener.UAV_U) * Lon_P * pi() / 180;
            }

            tf2::Quaternion myQuaternion;
            myQuaternion.setRPY(-UAV_roll_cmd, -UAV_pitch_cmd, UAV_yaw_cmd);
            r.orientation.x = myQuaternion.getX();
            r.orientation.y = myQuaternion.getY();
            r.orientation.z = myQuaternion.getZ();
            r.orientation.w = myQuaternion.getW();
            r.thrust = 0.3;
            pub.publish(r);
            modechange = false;
        } else {
            
            s.request.param_id = "SERVO3_MIN";
            s.request.value.integer = 1100;
            e.request.param_id = "SERVO4_FUNCTION";
            e.request.value.integer = 21;
            modechange = true;
           
        }
        pubth.call(s);
        pubth.call(e);

        ros::spinOnce();
        rate.sleep();
    
    }
    return 0;
}
