/*
 * Copyright (C) 2018 Yushu Yu and Yue Kang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "cluon-complete.hpp"
#include "opendlv-standard-message-set.hpp"

#include <cstdint>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <thread>
#include "dynamics_vehicle.h"
#include <math.h>


int32_t main(int32_t argc, char **argv) { 
    dynamics m_dynamics;    
    int32_t retCode{0};
    auto commandlineArguments = cluon::getCommandlineArguments(argc, argv);
    if ((0 == commandlineArguments.count("cid")) || (0 == commandlineArguments.count("freq")) || (0 == commandlineArguments.count("cid2") )) {
        std::cerr << argv[0] << " simulates can-rhino." << std::endl;
        std::cerr << "Usage:   " << argv[0] << " --cid=<OpenDaVINCI session> --cid2=<Second OD4Session> --freq=<Frequency> [--id=<Identifier in case of simulated units>] [--verbose]" << std::endl;
        std::cerr << "Example: " << argv[0] << " --cid=113 --cid2=114 --freq=50 [--save_file=/tmp/data_msg_nom_state.txt]" << std::endl;
        std::cerr << "(The second OD4Session is for vehicle state overriding from external source.)" << std::endl;
        retCode = 1;
    } else {
        const uint32_t ID{(commandlineArguments["id"].size() != 0) ? static_cast<uint32_t>(std::stoi(commandlineArguments["id"])) : 0};
        const bool VERBOSE{commandlineArguments.count("verbose") != 0};
        const int16_t FREQ = static_cast<uint16_t>(std::stoi(commandlineArguments["freq"]));
        const bool ifSave{commandlineArguments.count("save_file") != 0};
        auto filename{ifSave ? commandlineArguments["save_file"] : ""};

        // Interface to a running OpenDaVINCI session (ignoring any incoming Envelopes).
//        cluon::OD4Session od4{static_cast<uint16_t>(std::stoi(commandlineArguments["cid"])),
//            [&m_dynamics](auto){}
//        };
        uint16_t argCID = static_cast<uint16_t>(std::stoi(commandlineArguments["cid"]));
        std::shared_ptr<cluon::OD4Session> od4 = std::shared_ptr<cluon::OD4Session>(new cluon::OD4Session(argCID));

        uint16_t argCID2 = static_cast<uint16_t>(std::stoi(commandlineArguments["cid2"]));
        std::shared_ptr<cluon::OD4Session> od4_2 = std::shared_ptr<cluon::OD4Session>(new cluon::OD4Session(argCID2));

        // Define data triggered lambda functions
        auto Input_Pedal{[&m_dynamics, &VERBOSE](cluon::data::Envelope &&env)
            {
                // std::cout << "Now in lambda func. of Input_Pedal..." << std::endl;
                opendlv::proxy::PedalPositionRequest ppr = cluon::extractMessage<opendlv::proxy::PedalPositionRequest>(std::move(env));
                if (ppr.position() > 0.0f)
                {
                    m_dynamics.SetAcceleratorPedalPosition((double)ppr.position() * 100);
                    m_dynamics.SetBrakePedalPosition(0.0);
                    if (VERBOSE) std::cout << "Acc pedal request received: " << ppr.position() << std::endl;
                } else {
                    m_dynamics.SetAcceleratorPedalPosition(0.0);
                    m_dynamics.SetBrakePedalPosition((double)ppr.position() * 100);
                    if (VERBOSE) std::cout << "Break pedal request received: " << ppr.position() << std::endl;
                }
            }
        };

        auto Input_Steer{[&m_dynamics, &VERBOSE](cluon::data::Envelope &&env)
            {
                opendlv::proxy::GroundSteeringRequest gsr = cluon::extractMessage<opendlv::proxy::GroundSteeringRequest>(std::move(env));
                if (VERBOSE) std::cout << "Steering request received: " << gsr.groundSteering() << " (percentage)" << std::endl;
                if (gsr.groundSteering() > 1)
                {
                   // std::cerr << "WARNING: Right steering limit reached. Using max angle instead." << std::endl;
                    m_dynamics.SetRoadWheelAngle(m_dynamics.MAX_STEERING);
                }
                else if (gsr.groundSteering() < -1)
                {
                   // std::cerr << "WARNING: Left steering limit reached. Using max angle instead." << std::endl;
                       
                    m_dynamics.SetRoadWheelAngle(m_dynamics.MAX_STEERING * (-1));
                }
                else
                {
                    double angle = (double)gsr.groundSteering() * m_dynamics.MAX_STEERING;
                    m_dynamics.SetRoadWheelAngle(angle);
                }
            }
        };

        if (od4->isRunning())
        {
            od4->dataTrigger(opendlv::proxy::PedalPositionRequest::ID(), Input_Pedal);
            od4->dataTrigger(opendlv::proxy::GroundSteeringRequest::ID(), Input_Steer);
//            std::cout << "Data triggered functions attributed." << std::endl;
        }

        auto Input_Override_States{[&m_dynamics, &VERBOSE](cluon::data::Envelope &&env)
            {
                internal::nomState ext_state = cluon::extractMessage<internal::nomState>(std::move(env));
                if (VERBOSE) std::cout << "External overriding states received: " << std::endl \
                    << "[" << ext_state.xp_dot() << ", " << ext_state.yp_dot() << ", " << ext_state.psi_dot() << ", " << ext_state.epsi() << ", " \
                    << ext_state.ey() << ", " << ext_state.s() << ", " << ext_state.steer() << ", " << ext_state.acc() << "]" << std::endl;
                m_dynamics.state_global.epsi = ext_state.epsi();
                m_dynamics.state_global.ey = ext_state.ey();
                m_dynamics.state_global.s = ext_state.s();
                m_dynamics.state_global.v_body[0] = ext_state.xp_dot();
                m_dynamics.state_global.v_body[1] = ext_state.yp_dot();
                m_dynamics.state_global.omega_body[2] = ext_state.psi_dot(); // yaw rate
            }
        };
        /*message internal.nomState [id = 3002] {
  double xp_dot [id = 1];
  double yp_dot [id = 2];
  double psi_dot [id =3 ];
  double epsi [id = 4];
  double ey [id = 5];
  double s [id = 6];
  double steer [id = 7];
  double acc [id = 8];
}*/

        if (od4_2->isRunning())
        {
            od4_2->dataTrigger(internal::nomState::ID(), Input_Override_States);
        }


        // Define time triggered lambda function
        auto Output{[od4, &m_dynamics, &FREQ, &VERBOSE, &ifSave, &filename]() -> bool
            {
                uint16_t inner_freq = 1000 / FREQ;
                m_dynamics.T_samp = 0.001;  //sampling time
                for(uint16_t i = 0; i < inner_freq; i++)
                {
                    m_dynamics.integrator(VERBOSE);
                }

                opendlv::sim::KinematicState kinematicMsg;
                kinematicMsg.vx((float)m_dynamics.GetLongitudinalVelocity());
                kinematicMsg.vy((float)m_dynamics.GetLateralVelocity());
                kinematicMsg.vz(0.0f);
                kinematicMsg.rollRate(0.0f);
                kinematicMsg.pitchRate(0.0f);
                kinematicMsg.yawRate((float)m_dynamics.GetYawVelocity());
                od4->send(kinematicMsg);

                if (VERBOSE) std::cout << "Current Kinematic states sent." << std::endl;

                // 2018 Dec. Update: broadcast nominal states as well
                internal::nomState nomStateMsg;
                /*message internal.nomState [id = 3002] {
                  double xp_dot [id = 1];
                  double yp_dot [id = 2];
                  double psi_dot [id =3 ];
                  double epsi [id = 4];
                  double ey [id = 5];
                  double s [id = 6];
                  double steer [id = 7];
                  double acc [id = 8];
                }*/
                nomStateMsg.xp_dot(m_dynamics.GetLongitudinalVelocity());
                nomStateMsg.yp_dot(m_dynamics.GetLateralVelocity());
                nomStateMsg.psi_dot(m_dynamics.GetYawVelocity());
                nomStateMsg.epsi(m_dynamics.state_global.epsi);
                nomStateMsg.ey(m_dynamics.state_global.ey);
                nomStateMsg.s(m_dynamics.state_global.s);
                nomStateMsg.steer(m_dynamics.GetRoadWheelAngle());
                nomStateMsg.acc(m_dynamics.GetAcceleratorPedalPosition());

                od4->send(nomStateMsg);
                if (VERBOSE) std::cout << "Current nominal states sent." << std::endl;

                // Data saving into txt file (if "save_file" indicated)
                // number of rows = length of time 
                // each row contains the following data, seperated by tab: 
                // time nomStateMsg(all attributes)
                if (ifSave)
                {
                    std::ofstream txt(filename, std::ios::out | std::ios::app);
                    if (txt.is_open())
                    {
                        txt << ((double)clock())/CLOCKS_PER_SEC << '\t'
                            << nomStateMsg.xp_dot() << '\t' << nomStateMsg.yp_dot() << '\t' 
                            << nomStateMsg.psi_dot() << '\t' << nomStateMsg.epsi() << '\t' 
                            << nomStateMsg.ey() << '\t' << nomStateMsg.s() << '\t' 
                            << nomStateMsg.steer() << '\t' << nomStateMsg.acc() << '\n';
                        txt.close();
                    }
                    else std::cerr << "WARNING: Unable to save data into the file <" << filename << ">." << std::endl;
                }

                return false;

            }// end of lambda function
        }; // end of Output definition

        using namespace std::literals::chrono_literals;
        while (od4->isRunning()) {
            od4->timeTrigger(FREQ, Output);
//            std::this_thread::sleep_for(0.1s); // Commented as it is no longer simply data driven
//            std::cout << "Running..." << std::endl;
        }
    } // end of else (arguments check)
    return retCode;
}


dynamics::dynamics():
state_global({{2,0,0},{0,0,0},{0,0},0,0,0, 0, 0, 0}),
diff_global({{0,0,0},{0,0,0},{0,0},0,0,0,0,0,0}),
input_global({0,0,0,0}),
PI (3.14159265),
cp ( 20),  //parameter of cornering stiffness
mu ( 0.9),   //friction coefficient
mass ( 9840), //mass, FH16
g ( 9.8),  //acc due to gravity
rou ( 1.225),  C_d ( 0.7),  A ( 10),//coefficient of air drag, FH16
theta_g ( 0), //slope
lf ( 1.68),  //FH16
lr ( 1.715),   //FH16
Izz ( 41340),  //FH16
Je ( 4),  //flywheel inertia, FH16
a(1.68),
b(1.715),
cf (3.4812e+05),
cr (3.5537e+05),
m (9840),
Iz (41340),
psi_dot_com ( 0),
Efactor ( 0.5), ////fh16
i_final ( 3.46),  //FH16
i_gear(11.73),
eta_tr ( 1),  //efficient, FH16
eta_fd ( 0.9), //efficient, FH16
r_gear ( 0), //the flag of backward
a_xupper ( 1.35), //fh16
a_xlower ( 1.1),  //fh16
T_bmax ( 15000),
kb ( 2.25),
T_samp ( 0.001),
T_global ( 0),
agear ( 6),
agear_diff ( 0),
MAX_STEERING( 0.785398), // PI/4, 90 degree
omega_e ( 0),
Te ( 0)
{
	///////////////////the parameters of the vehicle//////////////////////

	PI = 3.14159265;
	//relate to wheels
	rw[0] = 0.435; rw[1] = 0.435;  //the radius of the wheel, FH16
	cp = 20;  //parameter of cornering stiffness
	mu = 0.9;   //friction coefficient
	Iw[0] = 11; Iw[1] = 11; //the inertia of the wheel, FH16
	fr[0] = 0.0164; fr[1] = 0.0164;  //rolling resistance coefficient, FH16

	mass = 9840; //mass, FH16
	g = 9.8;  //acc due to gravity
	rou = 1.225;  C_d = 0.7;  A = 10;//coefficient of air drag, FH16

	theta_g = 0; //slope
	lf = 1.68;  //FH16
	lr = 1.715;   //FH16
	Izz = 41340;  //FH16

	Je = 4;  //flywheel inertia, FH16

	//the fraction by which the engine torque is reduced
	Efactor = 0.5; ////fh16

	i_final = 3.46;  //FH16

	//fh16:
	i_tm[0] = 11.73; //gear ration for each gear
	i_tm[1] = 9.20; //gear ration for each gear
	i_tm[2] = 7.09; //gear ration for each gear
	i_tm[3] = 5.57; //gear ration for each gear
	i_tm[4] = 4.35; //gear ration for each gear
	i_tm[5] = 3.41; //gear ration for each gear
	i_tm[6] = 2.7; //gear ration for each gear
	i_tm[7] = 2.12; //gear ration for each gear
	i_tm[8] = 1.63; //gear ration for each gear
	i_tm[9] = 1.28; //gear ration for each gear
	i_tm[10] = 1.00; //gear ration for each gear
	i_tm[11] = 0.79; //gear ration for each gear


	//XC90
//	i_tm[0] = 5.2; //gear ration for each gear
//	i_tm[1] = 3.029; //gear ration for each gear
//	i_tm[2] = 1.96; //gear ration for each gear
//	i_tm[3] = 1.469; //gear ration for each gear
//	i_tm[4] = 1.231; //gear ration for each gear
//	i_tm[5] = 1; //gear ration for each gear
//	i_tm[6] = 0.809; //gear ration for each gear
//	i_tm[7] = 0.673; //gear ration for each gear
	//5.2,3.029,1.96,1.469,1.231,1,0.809,0.673
	// velocity_bound[2][10];

	i_gear = i_tm[0];
	eta_tr = 1;  //efficient, FH16
	eta_fd = 0.9; //efficient, FH16
	r_gear = 0; //the flag of backward

	//upper and lower bounds of acceleration limits,
	a_xupper = 1.35; //fh16
	a_xlower = 1.1;  //fh16

	//relating to brake:
	T_bmax = 15000;
	kb = 2.25;

	T_prop[0] = 0;  T_prop[1] = 0;
	T_brk[0] = 0;   T_brk[1] = 0;

	/////////////////////////////input//////////////////////////////
	input_global.A_ped = 0;   //pedal
	input_global.B_ped = 0;  //brake
	input_global.steering_angle = 0;

	/////////////////////////states/////////////////////////////////////
	//the velocity of the vehicle, expressed in the body frame of the vehicle
	for (int i = 0; i < 3; i++){
		state_global.v_body[i] = ((0 == i) ? 2 : 0);
		state_global.omega_body[i] = 0;  //angular velocity, body frame

		diff_global.vb_dot[i] = 0;  //dot of velocity of body
		diff_global.omegab_dot[i] = 0;  //dot of angular velocity of body
	}

	for (int i = 0; i < 2; i++){
	//the angular velocity of the wheel
		state_global.omega_w[i] = 0;
	//the derivative of angular velocity of the wheel
		diff_global.omega_wheel_dot[i] = 0;

	}
	//braking torque:
	state_global.T_b_general = 0;
	diff_global.T_b_dot_general = 0; //dot of T_b
	state_global.T_new_req = 0;
	state_global.Ttop = 0;

	//time step:
	T_samp = 0.001;
	T_global = 0;

	agear = 6;
	agear_diff = 0;

    MAX_STEERING = PI / 4;

	omega_e = 0;
 	Te = 0;

}



void dynamics::diff_equation(state_vehicle &state, input_vehicle &input,  double t_sim, diff_vehicle &out){
	//the differential equation of all the dynamics

	(void) t_sim;
        //int i = 0;
	//input:

        //double A_ped = input.A_ped;   //pedal
        //double B_ped = input.B_ped;  //brake
	double steering_angle = input.steering_angle;
        double acc_x = input.acc_x;

	/////////////////////////states/////////////////////////////////////
	//the velocity of the vehicle, expressed in the body frame of the vehicle
        //double v_body[3];
        //double omega_body[3];  //angular velocity, body frame
	//the angular velocity of the wheel
        //double omega_w[2];
	//braking torque:
        //double T_b_general;
        //double Ttop;
        //double T_new_req;
//	double vb_dot[3];  //dot of velocity of body
//	//the derivative of angular velocity of the wheel
//	double omegab_dot[3];  //dot of angular velocity of body
//	double omega_wheel_dot[2];
//	double T_b_dot_general; //dot of T_b


        /*
	//state:
	for (i = 0; i < 3; i++){
		v_body[i] = state.v_body[i];
		omega_body[i] = state.omega_body[i];  //angular velocity, body frame
	}
	for (i = 0; i < 2; i++){
			//the angular velocity of the wheel
		omega_w[i] = state.omega_w[i];
	}
			//braking torque:
		T_b_general = state.T_b_general;
		Ttop = state.Ttop;
		T_new_req = state.T_new_req;

	//wheel
	double f_sxy[2];  //3.14, middle variable
	double Fxy[2];  //3.15, middle variable
	double Fz[2];  //force along the z direction
	double Fw[3][2];  //force expressed in wheel frame
	double Fv[3][2];  //force expressed in body frame
	double T_roll[2];
	//the velocity of the vehicle, expressed in the wheel frame
	double vb_wheel[3][2];
	//slip
	double sx[2];
	double sy[2];
	double sxy[2];

	double delta[2];  //the steering angle
	delta[0] = steering_angle;
	delta[1] = 0;
	for (i=0; i<2; i++){
		vb_wheel[0][i] = v_body[0]*cos(delta[i]) + v_body[1] * sin(delta[i]);
		vb_wheel[1][i] = -v_body[0]*sin(delta[i]) + v_body[1] * cos(delta[i]);

		sx[i] = -(vb_wheel[0][i] - rw[i] * omega_w[i] ) / max_dynamics(abs_dynamics(rw[i]*omega_w[i]), 0.01);
		sy[i] = -(vb_wheel[1][i] ) / max_dynamics(abs_dynamics(rw[i]*omega_w[i]), 0.01);
		sxy[i] = sqrt(sx[i]*sx[i] + sy[i]*sy[i]);

		f_sxy[i] = 2/PI*atan(2*cp*sxy[i]/PI);
		//Fz[i] = mass*g*cos(theta_g)/2;
		Fz[0] = mass*g*cos(theta_g)*lf/(lf+lr);
		Fz[1] = mass*g*cos(theta_g)*lr/(lf+lr);
		Fxy[i] = mu*Fz[i]*f_sxy[i];
//		Fw[0][i] =(rw[i]*omega_w[i] - vb_wheel[0][i])*Fxy[i] / sqrt((rw[i]*omega_w[i] - vb_wheel[0][i])*(rw[i]*omega_w[i]
//		             - vb_wheel[0][i]) + vb_wheel[1][i] * vb_wheel[1][i]);
//		Fw[1][i] = - vb_wheel[1][i]*Fxy[i] / sqrt((rw[i]*omega_w[i] - vb_wheel[0][i])*(rw[i]*omega_w[i] - vb_wheel[0][i]) + vb_wheel[1][i] * vb_wheel[1][i]);
		Fw[0][i] = Fxy[i]*sx[i]/max_dynamics(sxy[i],0.1);
		Fw[1][i] = Fxy[i]*sy[i]/max_dynamics(sxy[i],0.1);

		if (omega_w[i] > 0.00000){
			T_roll[i] = fr[i]*Fz[i]*rw[i];
		}
		else if (omega_w[i] < -0.00000){
			T_roll[i] = -fr[i]*Fz[i]*rw[i];
		}
		else
			T_roll[i] = 0;
		//T_roll[i] = 0;  //test

		//force actuated on body, body frame
		Fv[0][i] = cos(delta[i]) * Fw[0][i] - sin(delta[i])*Fw[1][i];
		Fv[1][i] = sin(delta[i]) * Fw[0][i] + cos(delta[i])*Fw[1][i];
	}

	//brake for XC90:
	double T_req = T_bmax*0.01*B_ped;

 	T_brk[0] = T_b_general;
 	T_brk[1] = T_b_general;

//	T_brk[0] = 0; //test
//	T_brk[1] = 0; //test
	for(int j = 0; j < 2; j++){
		if(omega_w[j] < 0)
			T_brk[j] = -abs_dynamics(T_brk[j]);
		else if(omega_w[j] > 0)
			T_brk[j] =  abs_dynamics(T_brk[j]);
	}

	//body:
	double F_d[2];
	double F_g[2];
	double Fx; //total force
	double Fy;

	//x direction, body frame
	if(v_body[0] > 0.00){
		F_d[0] = 0.5*rou*C_d*A*v_body[0]*v_body[0];
		F_g[0] = mass*g*sin(theta_g);
	}
	else if(v_body[0] < -0.00){
		F_d[0] = -0.5*rou*C_d*A*v_body[0]*v_body[0];
		F_g[0] = mass*g*sin(theta_g);
	}
	else{
		F_d[0] = 0;
		F_g[0] = mass*g*sin(theta_g);
	}
	Fx = Fv[0][0] + Fv[0][1] - F_d[0] -F_g[0];

	//Fx = Fv[0][0] + Fv[0][1]; //test


	//y direction, body frame
//	F_d[1] = 0.5*rou*C_d*A*v_body[1]*v_body[1];
//	F_g[1] = mass*g*sin(theta_g);
//	Fy = Fv[1][0]+Fv[1][1] - F_d[1] -F_g[1];
	Fy = Fv[1][0] + Fv[1][1];

	double ax, ay;
	ax = Fx/mass;
	ay = Fy/mass;

	int drive_flag = 1;  // 1: rear driving, 0: rear and front driving
	//power strain, XC90:
	double omega_d, omega_d_dot;

	if(drive_flag == 0)
	{
		omega_d =  (omega_w[0] + omega_w[1])/2;
		omega_d_dot = (diff_global.omega_wheel_dot[0] + diff_global.omega_wheel_dot[1])/2;
	}

	if(drive_flag == 1){
		omega_d =  omega_w[1];
		omega_d_dot =   diff_global.omega_wheel_dot[1];
	}

	double omega_f, omega_f_dot;
	omega_f = omega_d*i_final;
	omega_f_dot = omega_d_dot*i_final;

	double omega_p, omega_p_dot;
	omega_p = omega_f;
	omega_p_dot = omega_f_dot;

	double omega_t, omega_t_dot;
	i_gear = i_tm[agear-1];  //agear is from 1 t0 ...,
	//i_gear = i_tm[0];  //test
	omega_t = omega_p * i_gear;
	omega_t_dot = omega_p_dot*i_gear;

	double omega_c, omega_c_dot;
	omega_c = omega_t;
	omega_c_dot = omega_t_dot;

	double  omega_e_dot;
	omega_e = omega_c;
	omega_e_dot =  omega_c_dot;

	double T_emax;

	//A_ped = 10; //test
  //  A_ped = T_global * 10;  //test
	if  ( (omega_e < 20) && (A_ped > 0))
	{
//		Teaped=10;  //test
//		T_emax = 10;  //test
		omega_e = 20;  //need initial speed to generate torque at initial time
	}

	T_emax = CalcEngineMaxTorque(omega_e);  //3.30

	double Teaped;
	Teaped = A_ped*0.01*T_emax;

	double T_alim;
	if(diff_global.vb_dot[0] > a_xupper)
		T_alim = Efactor*Teaped;
	else if (diff_global.vb_dot[0] < a_xlower)
		T_alim = Teaped;
	else
		T_alim = Teaped*((diff_global.vb_dot[0]-a_xlower)*(Efactor-1)/(a_xupper - a_xlower)+1);

	out.T_new_req_dot = Teaped - T_new_req;

	double T_req_alim = min_dynamics(T_new_req, T_alim);
	double Tsplit = 700;
	double Tbase = min_dynamics(T_req_alim, Tsplit);
	double Tdynreq = T_req_alim - Tbase;
	double k = 0.5;
	out.Ttop_dot = k*(Tdynreq - Ttop);
	Te = Tbase + Ttop;

	//(3.30)-(3.32):
//	if(diff_global.vb_dot[0] > a_xupper)
//		Te = Efactor*Teaped;
//	else if (diff_global.vb_dot[0] < a_xlower)
//		Te = T_emax;
//	else
//		Te = Teaped*((diff_global.vb_dot[0]-a_xlower)*(Efactor-1)/(a_xupper - a_xlower)+1);

	//T_emax = 700;
	//Te = T_emax;

	//omega_e_dot = (Te-Tc)/Je;
	double Tc;
	Tc = Te - Je*abs_dynamics(omega_e_dot);
	Tc = Te;

	//Tc = Te - Je*(omega_e_dot);

//	Te = 100;  //test
//	Tc = Te; //test
//	Tc = Te - Je*omega_e_dot;

	double k_speed_wtoe = i_final*i_gear;
	double	k_tau_ctow = eta_tr*i_gear*eta_fd*i_final;
//	double f2[2];
//	for (i=0; i<2; i++){
//    	f2[i] =  Fw[0][i] * rw[i] + T_roll[i] + T_brk[i];
//	}
	//Tc = (Iw[1] * Te + k_speed_wtoe*Je*f2[1])/(Iw[1]  + k_speed_wtoe*Je*k_tau_ctow); //according to Je*\dot omega_e = Te - Tc
	//Tc = Te;

	double Tt = Tc;
	double Tp = eta_tr*Tt*i_gear;
	double Tf = Tp;
	double Td = Tf*eta_fd*i_final;

	//XC90:
	double Twf;
	double Twr;

	//int drive_flag = 1;  // 1: rear driving, 0: rear and front driving
	if(drive_flag == 0){
		Twf = 0.4*Td;
		Twr = 0.6*Td;
	}
	if(drive_flag == 1){
		Twf = 0;
		Twr = 1*Td;
	}

	int rc = (r_gear > 0);

	double rev_trq = -1*rc*max_dynamics(Tt,0)*i_gear*i_final;

	if (r_gear == 1){
		T_prop[0] = 0.4*rev_trq;
		T_prop[1] = 0.6*rev_trq;
	}
	else{
		T_prop[0] = max_dynamics(Twf, 0);
		T_prop[1] = max_dynamics(Twr, 0);
	}

	//derivative part:
	double Te_direct[2];
	Te_direct[0] = 0;
	Te_direct[1] = Te;

//	Te_direct[0] = Te;
//	Te_direct[1] = Te;

	//derivative of wheel rotational velocity
	for (i=0; i<2; i++){
		double forceinducedtorque;
		forceinducedtorque = Fw[0][i] * rw[i];

		//out.omega_wheel_dot[i]= (T_prop[i] - T_brk[i] - forceinducedtorque - T_roll[i])/Iw[i];

     	//Je = 0; //if the flywheel is not considered
		out.omega_wheel_dot[i]=  (Te_direct[i] - 1/(k_tau_ctow)*(T_brk[i] + forceinducedtorque + T_roll[i] ))
				/(Je*k_speed_wtoe + Iw[i]/k_tau_ctow);

	}

	//derivative of body rotational velocity
	out.omegab_dot[2] = (lf*Fv[1][0] - lr* Fv[1][1])/Izz;
	//derivative of body velocity, expressed in body frame:
	out.vb_dot[0] = ax + v_body[1]*omega_body[2];
	out.vb_dot[1] = ay - v_body[0]*omega_body[2];
    //dot of T_b:
	out.T_b_dot_general = kb*(T_req-T_b_general);

	//agear
	double velocity_bound[2][12] =
	{{2.31,3,3.81,4.05,4.25,5.05,8.25,9.85,13.45,15.78,21.1,1000000},
	 {-1000000,2.2,2.28,3.5,3.9,4.1,4.77,6.1,6.5,8.5, 14.5,20}};

	if(mass > 12000){

		if ((v_body[0] > velocity_bound[0][agear-1]) && (agear <=11) )//upper bound
			agear_diff = 1;
		else if ( (v_body[0] < velocity_bound[1][agear-1]) && (agear >= 2))  //lower bound
			agear_diff = -1;
		else
			agear_diff = 0;
	}
	else{
		if ((v_body[0] > velocity_bound[0][agear-1]) && (agear <=11) )//upper bound
			agear_diff = 1;
		else if ( (v_body[0] < velocity_bound[1][agear-1]) && (agear >= 7))  //lower bound
			agear_diff = -1;
		else
			agear_diff = 0;
	}
        */

        //dynamics in the paper, input is steer angle and x-acc
        double xp_dot = state.v_body[0];
        double yp_dot = state.v_body[1];
        double psi_dot = state.omega_body[2];
        double epsi = state.epsi;
       // double ey = state.ey;
        // double s = state.s;


        if (state.v_body[0]> 1e-1)
        {
            //linerized model:
            out.vb_dot[0] = yp_dot*psi_dot + acc_x;   //dot xp_dot
            out.vb_dot[1] = -xp_dot*psi_dot -2*(cf+cr)/(m*xp_dot)*yp_dot-2*(a*cf-b*cr)/m/xp_dot*psi_dot + 2*cf/m*steering_angle;   // dot yp_dot
            out.omegab_dot[2] = -2*(a*cf-b*cr)/Iz/xp_dot*yp_dot-2*(a*a*cf+b*b*cr)/Iz/xp_dot*psi_dot + 2*a*cf/Iz*steering_angle;   //dot  psi_dot
            out.epsi_dot =  psi_dot - psi_dot_com;    // dot epsi
            out.ey_dot =  yp_dot*cos(epsi) + xp_dot*sin(epsi);    // dot ey
            out.s_dot =  xp_dot*cos(epsi)-yp_dot*sin(epsi);     // dot s
        }
        else
        {
            //linerized model:
            out.vb_dot[0] = 0;   //dot xp_dot
            out.vb_dot[1] =  0;   // dot yp_dot
            out.omegab_dot[2] =  0;   //dot  psi_dot
            out.epsi_dot =  0;    // dot epsi
            out.ey_dot =  0;    // dot ey
            out.s_dot =  0 ;     // dot s
        }



	////////test only, Feb. 12///
//    //parameters:
//    rw[0] = 0.347;
//    cp = 20;
//    mass = 2194;
//    g = 9.8;
//    theta_g=0;
//    mu=0.9;
//    double i_wheel = 11;
//    fr[0]= 0.0164;
//
//
//    sx[0] = -(v_body[0] - rw[0]  * omega_w[0] ) / max_dynamics(abs_dynamics(rw[0] *omega_w[0] ), 0.01);
//    sxy[0] = abs_dynamics(sx[0]);
//    f_sxy[0] = 2/PI*atan(2*cp*sxy[0]/PI);
//    double Fzz = mass*g*cos(theta_g)/2;
//    double Ftest = mu*Fzz *f_sxy[0];
//    Fw[0][0] = Ftest*sx[0]/max_dynamics(sxy[0],0.1);
//
//    out.vb_dot[0] = Fw[0][0]/mass;
//    out.omega_wheel_dot[0] = (50-Fw[0][0]*rw[0] )/i_wheel;

    /////////////////test finish/////////////////



	//ROS_INFO_STREAM("received path commands, flag_pc_cmd is set to)"<<vb_dot[0]));


//	std::cerr << "omega_w[0]: " << omega_w[0]  << "  omega_w[1]: " << omega_w[1] << std::endl;
//
//	std::cerr << "v_body[0]: " << v_body[0]
//	                                     << "	v_body[1]: " << v_body[1]
//<< "	omega_body[2]: " << omega_body[2] << std::endl;


//	std::cerr << "omega_wheel_dot[0]: " << out.omega_wheel_dot[0] <<
//			"  omega_wheel_dot[1]: " << out.omega_wheel_dot[1] << std::endl;
//	std::cerr << "vb_dot[0]: " << out.vb_dot[0] <<
//			"  vb_dot[1]: " << out.vb_dot[1] << std::endl;
//	std::cerr << "omegab_dot[2]: " << omegab_dot[2] << std::endl;
////
//
//	std::cerr << "Fw[0][0]: " << Fw[0][0]
//	          << "  Fw[1][0]: " << Fw[1][0]
//      << "   Fw[0][1]: " << Fw[0][1]<<
//      "   Fw[1][1]: " << Fw[1][1] << std::endl;
//
//	std::cerr << "sx[0]:" << sx[0]
//			<< "	sx[1]:" << sx[1]
//         << "	sy[0]:" << sy[0]
//           << "	sy[1]:" << sy[1] << std::endl;
//
//	std::cerr  << "Te: " << Te << "	Temx: " << T_emax << "  Engine speed: " << omega_e
//			<<"  speed f: " << omega_f<<std::endl;
//
//	std::cerr << "T_emax:" << T_emax << std::endl;
//
//	std::cerr << "T_prop[0]: " << T_prop[0]
// << "	T_prop[1]: " << T_prop[1]
//<< "	T_brk[0]: " << T_brk[0]
// << "	T_brk[1]: " << T_brk[1]
//<< "	T_roll[0]: " << T_roll[0]
//<< "	T_roll[1]: " << T_roll[1] << std::endl;
//


}



void dynamics::integrator(bool verbose){

	//update stete:
	int flag = 0;
	switch (flag)
	{
		{
			//4-order:
		case 0:
			state_vehicle  x2, x3, x4;
			diff_vehicle k1, k2, k3, k4;

			diff_equation(state_global, input_global,  T_global, k1);

			x2.T_b_general = state_global.T_b_general+ 0.5*T_samp*k1.T_b_dot_general;
			x2.Ttop = state_global.Ttop+ 0.5*T_samp*k1.Ttop_dot;
			x2.T_new_req = state_global.T_new_req+ 0.5*T_samp*k1.T_new_req_dot;
                        x2.epsi = state_global.epsi + 0.5*T_samp*k1.epsi_dot;
                        x2.ey = state_global.ey + 0.5*T_samp*k1.ey_dot;
                        x2.s = state_global.s + 0.5*T_samp*k1.s_dot;
			for(int i = 0; i < 3; i++){
				x2.v_body[i] = state_global.v_body[i] + 0.5*k1.vb_dot[i]*T_samp;
				x2.omega_body[i] = state_global.omega_body[i] + 0.5*k1.omegab_dot[i]*T_samp;
			}
			for(int j = 0; j < 2; j++){
				x2.omega_w[j] = 0.5*k1.omega_wheel_dot[j]*T_samp + state_global.omega_w[j];
			}
			diff_equation(x2, input_global,  T_global+0.5*T_samp, k2);

			x3.T_b_general = state_global.T_b_general+ 0.5*T_samp*k2.T_b_dot_general;
			x3.Ttop = state_global.Ttop+ 0.5*T_samp*k2.Ttop_dot;
			x3.T_new_req = state_global.T_new_req+ 0.5*T_samp*k2.T_new_req_dot;
                        x3.epsi = state_global.epsi + 0.5*T_samp*k2.epsi_dot;
                        x3.ey = state_global.ey + 0.5*T_samp*k2.ey_dot;
                        x3.s = state_global.s + 0.5*T_samp*k2.s_dot;
			for(int i = 0; i < 3; i++){
				x3.v_body[i] = state_global.v_body[i] + 0.5*k2.vb_dot[i]*T_samp;
				x3.omega_body[i] = state_global.omega_body[i] + 0.5*k2.omegab_dot[i]*T_samp;
			}
			for(int j = 0; j < 2; j++){
				x3.omega_w[j] = 0.5*k2.omega_wheel_dot[j]*T_samp + state_global.omega_w[j];
			}
			diff_equation(x3, input_global,  T_global + 0.5*T_samp, k3);

			x4.T_b_general = state_global.T_b_general+ T_samp*k3.T_b_dot_general;
			x4.Ttop = state_global.Ttop+ T_samp*k3.Ttop_dot;
			x4.T_new_req = state_global.T_new_req+  T_samp*k3.T_new_req_dot;
                        x4.epsi = state_global.epsi + 0.5*T_samp*k3.epsi_dot;
                        x4.ey = state_global.ey + 0.5*T_samp*k3.ey_dot;
                        x4.s = state_global.s + 0.5*T_samp*k3.s_dot;
			for(int i = 0; i < 3; i++){
				x4.v_body[i] = state_global.v_body[i] + k3.vb_dot[i]*T_samp;
				x4.omega_body[i] = state_global.omega_body[i] + k3.omegab_dot[i]*T_samp;
			}
			for(int j = 0; j < 2; j++){
				x4.omega_w[j] = k3.omega_wheel_dot[j]*T_samp + state_global.omega_w[j];
			}
			diff_equation(x4, input_global,  T_global + T_samp, k4);


			state_global.T_b_general = state_global.T_b_general+ T_samp*(k1.T_b_dot_general +
			2*k2.T_b_dot_general + 2*k3.T_b_dot_general+ k4.T_b_dot_general)/6;
			state_global.Ttop = state_global.Ttop+ T_samp*(k1.Ttop_dot +
					2*k2.Ttop_dot + 2*k3.Ttop_dot+ k4.Ttop_dot)/6;
			state_global.T_new_req = state_global.T_new_req+  T_samp*(k1.T_new_req_dot +
					2*k2.T_new_req_dot + 2*k3.T_new_req_dot+ k4.T_new_req_dot)/6;
			for(int i = 0; i < 3; i++){
				state_global.v_body[i] = state_global.v_body[i] + T_samp*(
						k1.vb_dot[i] + 2*k2.vb_dot[i] + 2*k3.vb_dot[i] + k4.vb_dot[i])/6;
				state_global.omega_body[i] = state_global.omega_body[i] + (
						k1.omegab_dot[i] + 2* k2.omegab_dot[i] + 2*k3.omegab_dot[i] + k4.omegab_dot[i])*T_samp/6;
			}
			for(int j = 0; j < 2; j++){
				state_global.omega_w[j] = state_global.omega_w[j] +  (
						k1.omega_wheel_dot[j] + 2*k2.omega_wheel_dot[j] + 2* k3.omega_wheel_dot[j] + k4.omega_wheel_dot[j])*T_samp/6;
			}
                        state_global.epsi = state_global.epsi+T_samp*(k1.epsi_dot + 2*k2.epsi_dot + 2*k3.epsi_dot+ k4.epsi_dot)/6;
                        state_global.ey = state_global.ey+T_samp*(k1.ey_dot + 2*k2.ey_dot + 2*k3.ey_dot+ k4.ey_dot)/6;
                        state_global.s = state_global.s+T_samp*(k1.s_dot + 2*k2.s_dot + 2*k3.s_dot+ k4.s_dot)/6;

			T_global = T_global+T_samp;
			diff_global.vb_dot[0] = ( k1.vb_dot[0] + 2*k2.vb_dot[0] + 2*k3.vb_dot[0] + k4.vb_dot[0])/6;
			diff_global.omegab_dot[0] = ( k1.omegab_dot[0] + 2*k2.omegab_dot[0] + 2*k3.omegab_dot[0] + k4.omegab_dot[0])/6;
			diff_global.omegab_dot[1] = ( k1.omegab_dot[1] + 2*k2.omegab_dot[1] + 2*k3.omegab_dot[1] + k4.omegab_dot[1])/6;
			diff_global.omega_wheel_dot[1] = ( k1.omega_wheel_dot[1] +
					2*k2.omega_wheel_dot[1] + 2*k3.omega_wheel_dot[1] + k4.omega_wheel_dot[1])/6;
                        diff_global.epsi_dot= ( k1.epsi_dot + 2*k2.epsi_dot + 2*k3.epsi_dot + k4.epsi_dot)/6;
                        diff_global.ey_dot= ( k1.ey_dot + 2*k2.ey_dot + 2*k3.ey_dot + k4.ey_dot)/6;
                        diff_global.s_dot= ( k1.s_dot + 2*k2.s_dot + 2*k3.s_dot + k4.s_dot)/6;

			break;
		}



		{ //1-order:
		case 1:

			diff_vehicle k_1st;
			diff_equation(state_global, input_global,  T_global, k_1st);
			state_global.T_b_general = state_global.T_b_general + T_samp*k_1st.T_b_dot_general;
			//body:
			for(int i = 0; i < 3; i++){
				state_global.v_body[i] = state_global.v_body[i] + k_1st.vb_dot[i]*T_samp;
				state_global.omega_body[i] = state_global.omega_body[i] + k_1st.omegab_dot[i]*T_samp;
			}
			//wheel:
			for(int j = 0; j < 2; j++){
				state_global.omega_w[j] = k_1st.omega_wheel_dot[j]*T_samp + state_global.omega_w[j];
			}
			T_global = T_global+T_samp;
			break;
		}

	}
	//agear
	agear = agear_diff + agear;

    if (verbose && false)
    {
        std::cerr << "angular velocity of wheel: " << state_global.omega_w[0]  << "," <<  state_global.omega_w[1] << std::endl;
    ////    std::cerr << "angular acc of wheel:  " << diff_global.omega_wheel_dot[0] << ","  << diff_global.omega_wheel_dot[1] << std::endl;
    //
        std::cerr << "body velocity (x, y, rot z): " << state_global.v_body[0] << ","
                << state_global.v_body[1] << "," << state_global.omega_body[2] << std::endl;
    ////    std::cerr << "body acc (x, y, rot z): " << diff_global.vb_dot[0]  << ","
    ////            << diff_global.vb_dot[1]  << "," << diff_global.omegab_dot[2] << std::endl;
    //
        std::cerr << "i_gear: " << i_gear
    << "    a_gear: " << agear << " diff_gear: " << agear_diff << std::endl;

        std::cerr << "Time: "   << T_global << std::endl;

        std::cerr << std::endl;
    }



}


double dynamics::max_dynamics(double aa, double bb){
        if (aa>=bb)
                return aa;
	else
                return bb;

}

double dynamics::min_dynamics(double aa, double bb){
        if (aa<=bb)
                return aa;
	else
                return bb;

}

double dynamics::abs_dynamics(double aa){
        if (aa>=0)
                return aa;
	else{
                double bb = -aa;
                return bb;
	}

}

double dynamics::CalcEngineMaxTorque(double m_engineSpeed) {
	int size = 17;
	double torqueLookupTable[17][2] =
	{	  {0.0, 0.0},
			  {62.8318, 1660.0},
			  {73.3038, 1880.0},
			  {83.775, 2240.0},
			  {94.247, 2900.0},
			  {104.719, 3550.0},
			  {115.191, 3550.0},
			  {125.663, 3550.0},
			  {136.135, 3550.0},
			  {146.607, 3550.0},
			  {157.079, 3470.0},
			  {167.551, 3310.0},
			  {178.023, 3120.0},
			  {188.495, 2880.0},
			  {198.967, 2660.0},
			  {209.439, 1680.0},
			  {219.911, 0.0}
	};

	if (m_engineSpeed < torqueLookupTable[0][0]) {
	return 0.0;
	}

	if (m_engineSpeed >= torqueLookupTable[size- 1][0]) {
	return 0.0;
	}

	for (int i = 0; i < size - 1; i++) {
		double const x1 = torqueLookupTable[i][0];
		double const x2 = torqueLookupTable[i+1][0];

		if (m_engineSpeed >= x1 && m_engineSpeed < x2) {
 	      double const r = (m_engineSpeed - x1) / (x2 - x1);

    	  double const y1 = torqueLookupTable[i][1];
		  double const y2 = torqueLookupTable[i+1][1];
 		  double const maxTorque = y1 + r * (y2 - y1);

		  return maxTorque;
		}
	}

	std::cerr << "Lookup failed. This should never happen." << m_engineSpeed <<  std::endl;
	return 0.0;
}

double dynamics::GetLongitudinalVelocity() const{
	return state_global.v_body[0];
}

double dynamics::GetLateralAcceleration() const{
	return diff_global.vb_dot[1];

}
double dynamics::GetLateralVelocity() const{
	return state_global.v_body[1];

}
double dynamics::GetLongitudinalAcceleration() const{
	return diff_global.vb_dot[0];
}

double dynamics::GetYawAcceleration() const{
	return diff_global.omegab_dot[2];
}
double dynamics::GetYawVelocity() const{
	return state_global.omega_body[2];
}

double dynamics::GetAcceleratorPedalPosition() const{
	return input_global.A_ped;
}

double dynamics::GetEngineSpeed() const{
	return omega_e;
}

double dynamics::GetEngineTorque() const{
	return Te;
}

void dynamics::SetAcceleratorPedalPosition(double pos){
	input_global.A_ped = pos;
}

double dynamics::GetBrakePedalPosition() const{
	return input_global.B_ped;
}

void dynamics::SetBrakePedalPosition(double pos){
	input_global.B_ped = pos;
}

int32_t dynamics::GetGear() const{
	return agear;
}

double dynamics::GetFrontWheelSpeed() const{
	return state_global.omega_w[0];
}
double dynamics::GetRearWheelSpeed() const{
	return state_global.omega_w[1];
}
double dynamics::GetRoadWheelAngle() const{
	return input_global.steering_angle;
}
void dynamics::SetRoadWheelAngle(double wa){
        input_global.steering_angle = wa;
}




