#include <stddef.h>
#include <stdio.h>    
#include <math.h> 
#include <cmath>
#include <sstream>             
#include <iostream>
#include <tf2/LinearMath/Quaternion.h>
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

using namespace std;
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }
double pix_Y2mpc_Y(double x) {
	 if(x-pi()/2<-pi())
	 {
		return -(x-pi()/2+2*pi());
	 }
	 else
	 {
		return -(x-pi()/2);
	 }
	 }

class Listener
{
	public:
    double qua_w, qua_x, qua_y, qua_z;
    double siny_cosp, cosy_cosp;
	double UAV_yaw;
	double UAV_E, UAV_N, UAV_U, UAV_p, UAV_q, UAV_r;
	double  guided_state;

	void imu_callback(const sensor_msgs::Imu::ConstPtr& msg);
	void pos_callback(const nav_msgs::Odometry::ConstPtr& msg);
	void guided_callback(const mavros_msgs::State::ConstPtr& msg);
};

void Listener::imu_callback(const sensor_msgs::Imu::ConstPtr& msg) 
{
  std_msgs::String pub_str;
  std::stringstream ss;
  ss << "controller heard: x: " << msg->orientation.x << " y: " << msg->orientation.y << " z: " << msg->orientation.z << " w: " << msg->orientation.w;
  pub_str.data = ss.str();

  qua_w = msg->orientation.w;
  qua_x = msg->orientation.x;
  qua_y = msg->orientation.y;
  qua_z = msg->orientation.z;

 UAV_p = msg->angular_velocity.x;
 UAV_q= msg->angular_velocity.y;
 UAV_r = msg->angular_velocity.z;

 siny_cosp = 2 * (qua_w * qua_z + qua_x * qua_y);
 cosy_cosp = 1 - 2 * (qua_y * qua_y + qua_z * qua_z);
 UAV_yaw = std::atan2(siny_cosp, cosy_cosp);
}

void Listener::pos_callback(const nav_msgs::Odometry::ConstPtr& msg) 
{
	UAV_E = msg->pose.pose.position.x;   // East(x)
	UAV_N = msg->pose.pose.position.y;   // North(y)
	UAV_U = msg->pose.pose.position.z;   // Up(z)
}

void Listener::guided_callback(const mavros_msgs::State::ConstPtr& msg)
{
	guided_state = msg->guided;
}

// Instance of model class
static thesis_mpc rtObj;           
//    
// Path Choose [0:Straight Line, 1:Circular Path]
double path_type = 0; 
//Guidanse Law Parameters
double path_psi = 0;
double path_N_start= 0;
double path_E_start = 0;
double over_circle_N = 0;
double over_circle_E = 0;
double path_N_mid;
double path_E_mid;
double rotation_heading;
double radius = 100;
double gostraight = 30;
double L1 = 60;
double V = 15;
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
double ahead_uav_E ;
double ahead_uav_N ;
double a_cmd;
double distance_error;
double delta_psi;
double UAV_roll_cmd;
double UAV_yaw_cmd;
//
//Atitude Hold Parameter
double Up_cmd; 
double Lon_P;
double UAV_pitch_cmd;
//
int_T main(int_T argc, char **argv)
{
	ros::init(argc, argv, "thesis_main");
	ros::NodeHandle n;
	Listener listener;
	ros::Subscriber sub_imu = n.subscribe<sensor_msgs::Imu> ("mavros/imu/data", 1, &Listener::imu_callback, &listener);
	ros::Subscriber sub_position = n.subscribe<nav_msgs::Odometry> ("mavros/global_position/local", 1, &Listener::pos_callback, &listener);
	ros::Subscriber sub_guided = n.subscribe<mavros_msgs::State> ("mavros/state", 1, &Listener::guided_callback, &listener);
	mavros_msgs::AttitudeTarget r;
	mavros_msgs::ParamSet s;
	mavros_msgs::ParamSet e;
	ros::Publisher pub = n.advertise<mavros_msgs::AttitudeTarget>("mavros/setpoint_raw/attitude", 1); 
	ros::ServiceClient pubth = n.serviceClient<mavros_msgs::ParamSet>("mavros/param/set" );
	ros::Rate rate(10.0);
	
    while (ros::ok())
	{
		if(listener.guided_state == 1)
		{
			s.request.param_id = "SERVO3_MIN";
			s.request.value.integer = 1700;

			e.request.param_id = "SERVO4_FUNCTION";
			e.request.value.integer =0;
    		r.header.stamp = ros::Time::now();
			if(path_type==0)
			{
				//Guidance Law (Straight Line)
				if(path_psi ==0)
				{
					path_psi = pix_Y2mpc_Y(listener.UAV_yaw);
					path_N_start = listener.UAV_N;
					path_E_start = listener.UAV_E;
				}
				path_N_end = path_N_start + V*cos(path_psi);
				path_E_end = path_E_start + V*sin(path_psi);
				line_a = path_N_end - path_N_start;
				line_b = path_E_start - path_E_end;
				line_c = path_E_end*path_N_start - path_E_start*path_N_end;
				distance_error = (line_a*listener.UAV_E+line_b*listener.UAV_N+line_c)/(sqrt(pow(line_a,2)+pow(line_b,2)));
				if(abs((line_a*listener.UAV_E+line_b*listener.UAV_N+line_c)/(sqrt(pow(line_a,2)+pow(line_b,2)))) <= L1)
				{
					if(line_b == 0)
					{
						Targetpoint_1[1] = -line_c/line_a;
						Targetpoint_1[0] = (-(-2*listener.UAV_N)+sqrt(pow(-2*listener.UAV_N,2)-4*(2*listener.UAV_E*line_c/line_a+pow(line_c/line_a,2)-pow(L1,2)+pow(listener.UAV_N,2)+pow(listener.UAV_E,2))))/2;
						Targetpoint_2[1] = -line_c/line_a;
						Targetpoint_2[0] = (-(-2*listener.UAV_N)-sqrt(pow(-2*listener.UAV_N,2)-4*(2*listener.UAV_E*line_c/line_a+pow(line_c/line_a,2)-pow(L1,2)+pow(listener.UAV_N,2)+pow(listener.UAV_E,2))))/2;
					}
					else
					{
						Targetpoint_1[1] = (-(2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N)+sqrt(pow((2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N),2)-4*(pow((-line_a/line_b),2)+1)*(pow((-line_c/line_b),2)-2*listener.UAV_N*(-line_c/line_b)-pow(L1,2)+pow(listener.UAV_E,2)+pow(listener.UAV_N,2))))/(2*(pow((-line_a/line_b),2)+1)); //East
						Targetpoint_1[0] = (-line_a/line_b)*Targetpoint_1[1]+(-line_c/line_b);                                                                //North
						Targetpoint_2[1] = (-(2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N)-sqrt(pow((2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N),2)-4*(pow((-line_a/line_b),2)+1)*(pow((-line_c/line_b),2)-2*listener.UAV_N*(-line_c/line_b)-pow(L1,2)+pow(listener.UAV_E,2)+pow(listener.UAV_N,2))))/(2*(pow((-line_a/line_b),2)+1)); //East
						Targetpoint_2[0] = (-line_a/line_b)*Targetpoint_2[1]+(-line_c/line_b);                                                               //North
					}
					ahead_uav_E = listener.UAV_E + V*cos(listener.UAV_yaw);
					ahead_uav_N = listener.UAV_N + V*sin(listener.UAV_yaw);
					choose_Tp1 = sqrt(pow(ahead_uav_E-Targetpoint_1[1],2)+pow(ahead_uav_N-Targetpoint_1[0],2));
					choose_Tp2 = sqrt(pow(ahead_uav_E-Targetpoint_2[1],2)+pow(ahead_uav_N-Targetpoint_2[0],2));
					if(abs(choose_Tp1)<abs(choose_Tp2))
					{
						a_cmd = 2*pow(V,2)/L1*sin(pix_Y2mpc_Y(std::atan2(Targetpoint_1[0]-listener.UAV_N,Targetpoint_1[1]-listener.UAV_E))-pix_Y2mpc_Y(listener.UAV_yaw));
					}
					else
					{
						a_cmd = 2*pow(V,2)/L1*sin(pix_Y2mpc_Y(std::atan2(Targetpoint_2[0]-listener.UAV_N,Targetpoint_2[1]-listener.UAV_E))-pix_Y2mpc_Y(listener.UAV_yaw));
					}
				}
				else
				{
					if((line_a*listener.UAV_E+line_b*listener.UAV_N+line_c)/(sqrt(pow(line_a,2)+pow(line_b,2))) > 0)
					{
						a_cmd = 2*pow(V,2)/L1*sin(pi());
					}
					else
					{
						a_cmd = 2*pow(V,2)/L1*sin(-pi());
					}
				}
				delta_psi = a_cmd/V;
			}
			else
			{
				//Guidance Law (Circular Path)
				if(path_psi ==0)
				{
					path_psi = listener.UAV_yaw;
					path_N_start = listener.UAV_N;
					path_E_start = listener.UAV_E;
				}
				path_N_mid = path_N_start +30*sin(path_psi);
				path_E_mid = path_E_start+30*cos(path_psi);
				if (path_psi - pi()/2<-pi())
				{
					rotation_heading = path_psi - pi()/2 + 2*pi();
				}
				else
				{
					rotation_heading = path_psi - pi()/2;
				}
				path_N_end = path_N_mid + 100*sin(rotation_heading);
				path_E_end = path_E_mid + 100*cos(rotation_heading);
				distance_error = radius-sqrt(pow((path_E_end-listener.UAV_E),2)+pow((path_N_end-listener.UAV_N),2));
				if(distance_error<=L1)
				{
					Inter_d = sqrt(pow(abs(path_E_end-listener.UAV_E),2)+pow(abs(path_N_end-listener.UAV_N),2));
					Inter_A = (pow(L1,2)-pow(radius,2)+pow(Inter_d,2))/(2*Inter_d);
					Inter_h = sqrt(pow(L1,2)-pow(Inter_A,2));
					Inter_E = listener.UAV_E+Inter_A*(path_E_end-listener.UAV_E)/Inter_d;
					Inter_N = listener.UAV_N+Inter_A*(path_N_end-listener.UAV_N)/Inter_d;
					Targetpoint_1[1] = Inter_E-Inter_h*(path_N_end-listener.UAV_N)/Inter_d;
					Targetpoint_1[0] = Inter_N+Inter_h*(path_E_end-listener.UAV_E)/Inter_d;
					Targetpoint_2[1] = Inter_E+Inter_h*(path_N_end-listener.UAV_N)/Inter_d;
					Targetpoint_2[0] = Inter_N-Inter_h*(path_E_end-listener.UAV_E)/Inter_d;
				}
				else
				{
					if(over_circle_E ==0)
					{
						over_circle_N = listener.UAV_N;
						over_circle_E = listener.UAV_E;
					}
					line_a = path_N_end - over_circle_N ;
					line_b = over_circle_E - path_E_end;
					line_c = path_E_end*over_circle_N  - over_circle_E*path_N_end;
					if(line_b == 0)
					{
						Targetpoint_1[1] = -line_c/line_a;
						Targetpoint_1[0] = (-(-2*listener.UAV_N)+sqrt(pow(-2*listener.UAV_N,2)-4*(2*listener.UAV_E*line_c/line_a+pow(line_c/line_a,2)-pow(L1,2)+pow(listener.UAV_N,2)+pow(listener.UAV_E,2))))/2;
						Targetpoint_2[1] = -line_c/line_a;
						Targetpoint_2[0] = (-(-2*listener.UAV_N)-sqrt(pow(-2*listener.UAV_N,2)-4*(2*listener.UAV_E*line_c/line_a+pow(line_c/line_a,2)-pow(L1,2)+pow(listener.UAV_N,2)+pow(listener.UAV_E,2))))/2;
					}
					else
					{
						Targetpoint_1[1] = (-(2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N)+sqrt(pow((2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N),2)-4*(pow((-line_a/line_b),2)+1)*(pow((-line_c/line_b),2)-2*listener.UAV_N*(-line_c/line_b)-pow(L1,2)+pow(listener.UAV_E,2)+pow(listener.UAV_N,2))))/(2*(pow((-line_a/line_b),2)+1)); //East
						Targetpoint_1[0] = (-line_a/line_b)*Targetpoint_1[1]+(-line_c/line_b);                                                                //North
						Targetpoint_2[1] = (-(2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N)-sqrt(pow((2*(-line_a/line_b)*(-line_c/line_b) -2*listener.UAV_E-2*(-line_a/line_b)*listener.UAV_N),2)-4*(pow((-line_a/line_b),2)+1)*(pow((-line_c/line_b),2)-2*listener.UAV_N*(-line_c/line_b)-pow(L1,2)+pow(listener.UAV_E,2)+pow(listener.UAV_N,2))))/(2*(pow((-line_a/line_b),2)+1)); //East
						Targetpoint_2[0] = (-line_a/line_b)*Targetpoint_2[1]+(-line_c/line_b);                                                               //North
					}
				}
				ahead_uav_E = listener.UAV_E + V*cos(listener.UAV_yaw);
				ahead_uav_N = listener.UAV_N + V*sin(listener.UAV_yaw);
				choose_Tp1 = sqrt(pow(ahead_uav_E-Targetpoint_1[1],2)+pow(ahead_uav_N-Targetpoint_1[0],2));
				choose_Tp2 = sqrt(pow(ahead_uav_E-Targetpoint_2[1],2)+pow(ahead_uav_N-Targetpoint_2[0],2));
				if(abs(choose_Tp1)<abs(choose_Tp2))
				{
					a_cmd = 2*pow(V,2)/L1*sin(pix_Y2mpc_Y(std::atan2(Targetpoint_1[0]-listener.UAV_N,Targetpoint_1[1]-listener.UAV_E))-pix_Y2mpc_Y(listener.UAV_yaw));
				}
				else
				{
					a_cmd = 2*pow(V,2)/L1*sin(pix_Y2mpc_Y(std::atan2(Targetpoint_2[0]-listener.UAV_N,Targetpoint_2[1]-listener.UAV_E))-pix_Y2mpc_Y(listener.UAV_yaw));
				}
				delta_psi = a_cmd/V;
			}
			//Path Following Position Controller
			rtObj.rtU.uav_r = -listener.UAV_r ;
			rtObj.rtU.r_cmd = delta_psi;
			// for (int i = 0; i < 20; i++)
			// {
				rtObj.step();
			// }
			UAV_roll_cmd = rtObj.rtY.Out1;
			UAV_yaw_cmd = 0;
			//
			// Height Control
			Up_cmd = 70;
			Lon_P = 0.6;
			if((Up_cmd-listener.UAV_U)*Lon_P>10)
			{
				UAV_pitch_cmd = 10*pi()/180;
			}
			else if ((Up_cmd-listener.UAV_U)*Lon_P<-10)
			{
				UAV_pitch_cmd = -10*pi()/180;
			}
			else
			{
				UAV_pitch_cmd = (Up_cmd-listener.UAV_U)*Lon_P*pi()/180;
			}
			
			//
			tf2::Quaternion myQuaternion;
			myQuaternion.setRPY( UAV_roll_cmd, -UAV_pitch_cmd, UAV_yaw_cmd);
			r.orientation.x = myQuaternion.getX(); 
        	r.orientation.y = myQuaternion.getY();
        	r.orientation.z = myQuaternion.getZ();
        	r.orientation.w = myQuaternion.getW();
			r.thrust = 0.5;
    		pub.publish(r);
			
			if(path_type==0)
			{
				printf("in_x: %f,in_y: %f,in_yaw: %f,roll_cmd:%f \n",path_E_start,path_N_start,rad2deg(path_psi),rad2deg(UAV_roll_cmd));
			}
			else
			{
				printf("Cir_x: %f,Cir_y: %f,in_yaw: %f \n",path_E_end,path_N_end,rad2deg(path_psi));
			}
			printf("UAV_r: %f, r_cmd: %f \n",-listener.UAV_r,delta_psi);
			printf("Distance: %f \n",distance_error);
			printf("phi_cmd: %f \n",rad2deg(UAV_roll_cmd));
			if(abs(choose_Tp1)<abs(choose_Tp2))
			{
				printf("Target_E: %f , Target_N: %f \n",Targetpoint_1[1],Targetpoint_1[0]);
			}
			else
			{
				printf("Target_E: %f , Target_N: %f \n",Targetpoint_2[1],Targetpoint_2[0]);
			}
			printf(" ******************************* \n");
		}
		else
		{
			s.request.param_id = "SERVO3_MIN";
			s.request.value.integer = 1100;
			e.request.param_id = "SERVO4_FUNCTION";
			e.request.value.integer = 21;
			printf("%f  \n",listener.guided_state);
			if(path_type==0)
			{
				printf("in_x: %f,in_y: %f,in_yaw: %f \n",path_E_start,path_N_start,rad2deg(path_psi));
			}
			else
			{
				printf("Cir_x: %f,Cir_y: %f,in_yaw: %f \n",path_E_end,path_N_end,rad2deg(path_psi));
			}
			printf(" ******************************* \n");
		}	
		pubth.call(s);
		pubth.call(e);


		//pub1.publish(s);
		ros::spinOnce();
		rate.sleep();
	}
	return 0;
}



