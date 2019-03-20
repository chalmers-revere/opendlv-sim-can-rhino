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

#include <cstdint>
#include <chrono>
#include <iostream>
#include <sstream>
#include <thread>

#include <cmath>
#include <Eigen/Dense>

#include "cluon-complete.hpp"
#include "opendlv-standard-message-set.hpp"

#include "data_structure.hpp"

int32_t main(int32_t argc, char *argv[])
{
    bool flag_nomu(true), flag_nomstate(true), flag_actstate(true); 
    auto commandlineArguments = cluon::getCommandlineArguments(argc, argv);

    uint16_t CID = 113, CID2 = 114;
    if (0 == commandlineArguments.count("cid"))
    {
        std::cerr << "WARNING: No cid assigned, using 113 by default ('--cid=113')." << std::endl;
    }
    else CID = std::stoi(commandlineArguments["cid"]);

    if (0 == commandlineArguments.count("cid2"))
    {
        std::cerr << "WARNING: No cid2 assigned, using 114 by default ('--cid2=114')." << std::endl;
    }
    else CID2 = std::stoi(commandlineArguments["cid2"]);

    uint32_t FREQ = 50;
    if (0 == commandlineArguments.count("freq"))
    {
        std::cerr << "WARNING: No frequency assigned, using 50 by default ('--freq=50')." << std::endl;
    }
    else FREQ = std::stoi(commandlineArguments["freq"]);

    double k_scale_steer = 1.0;
    if (0 == commandlineArguments.count("k_scale_steer"))
    {
        std::cerr << "WARNING: No k_scale_steer assigned, using 50 by default ('--k_scale_steer=1.0')." << std::endl;
    }
    else k_scale_steer = std::stoi(commandlineArguments["k_scale_steer"]);
    double k_scale_acc = 1.0;
    if (0 == commandlineArguments.count("k_scale_acc"))
    {
        std::cerr << "WARNING: No k_scale_acc assigned, using 50 by default ('--k_scale_acc=1.0')." << std::endl;
    }
    else k_scale_acc = std::stoi(commandlineArguments["k_scale_acc"]);


    bool const VERBOSE{commandlineArguments.count("verbose") != 0};



    FB_state nom_state(16.0, 0, 0, 0, 0, 0, 0, 0);
    if (VERBOSE) 
        std::cout << "nom_state initialised." << std::endl;

    FB_state real_state(16.0, 0, 0, 0, 0, 0, 0, 0);

    Eigen::Vector2d nom_u; // Order: [acc, steer]'
    nom_u << 0.0, 0.0;
    if (VERBOSE) 
        std::cout << "nom_u initialised." << std::endl;

    // Global_variables gl;

    cluon::OD4Session od4(CID, [&nom_state, &nom_u, &VERBOSE, &flag_nomu, &flag_nomstate](cluon::data::Envelope &&env) noexcept {
        if (env.dataType() == internal::nomState::ID())
        {
            internal::nomState received = cluon::extractMessage<internal::nomState>(std::move(env));
            nom_state = FB_state(received.xp_dot(), received.yp_dot(), received.psi_dot(), received.epsi(), received.ey(), received.s(), received.steer(), received.acc());
            flag_nomstate = false; 
            if (VERBOSE)
            {
                std::cout << "New nom_state received:" << std::endl;
                nom_state.print(); 
            }
        }
        else if (env.dataType() == internal::nomU::ID())
        {
            internal::nomU received = cluon::extractMessage<internal::nomU>(std::move(env));
            nom_u << received.steer(), received.acc();  //notice the order
            flag_nomu = false; 
            if (VERBOSE)
            {
                std::cout << "New nom_u received:" << std::endl;
                std::cout << "[" << nom_u(0) << ", " << nom_u(1) << "]" << std::endl;
            }
        }
    });
    cluon::OD4Session od4_2(CID2, [&real_state, &VERBOSE, &flag_actstate](cluon::data::Envelope &&env) noexcept {
        if (env.dataType() == opendlv::sim::Frame::ID())
        {
            opendlv::sim::Frame received = cluon::extractMessage<opendlv::sim::Frame>(std::move(env));
            real_state.s = received.x();
            real_state.ey = received.y();
            real_state.epsi = received.yaw();

            if (VERBOSE)
            {
                std::cout << "Vehicle position updated:" << std::endl;
                std::cout << "x:" << real_state.s << ", y:" << real_state.ey << ", yaw:" << real_state.epsi << std::endl;
            }
        }
        else if (env.dataType() == opendlv::sim::KinematicState::ID())
        {
            opendlv::sim::KinematicState received = cluon::extractMessage<opendlv::sim::KinematicState>(std::move(env));
            real_state.xp_dot = received.vx();
            real_state.yp_dot = received.vy();
            real_state.psi_dot = received.yawRate();
            if (VERBOSE)
            {
                std::cout << "Vehicle velocity updated:" << std::endl;
                std::cout << "vx:" << real_state.xp_dot << ", y:" << real_state.yp_dot << ", yawRate:" << real_state.psi_dot << std::endl;
            }
        }
        flag_actstate = false; 
        /*if (env.dataType() == internal::nomState::ID())
        {
            internal::nomState received = cluon::extractMessage<internal::nomState>(std::move(env));
            real_state = FB_state(received.xp_dot(), received.yp_dot(), received.psi_dot(), received.epsi(), received.ey(), received.s(), received.steer(), received.acc());
            if (VERBOSE)
            {
                std::cout << "New actual state received:" << std::endl;
                real_state.print(); 
            }
        }*/
    });
    if (0 == od4.isRunning())
    {
        std::cerr << "ERROR: Internal OD4 not running!!!" << std::endl;
        return -1;
    }
    if (0 == od4_2.isRunning())
    {
        std::cerr << "ERROR: External OD4 not running!!!" << std::endl;
        return -2;
    }
    while (od4.isRunning() && od4_2.isRunning())
    {
        auto sendMsg{[&od4_2, &nom_state, &real_state, &nom_u, &VERBOSE, &flag_nomu, &flag_nomstate, &flag_actstate, 
                      &k_scale_steer, &k_scale_acc]() -> bool
            {
                Eigen::Vector2d u;
                if (nom_state.xp_dot <= 1e-1)
                {
                    u << 0, -5;
                    if (VERBOSE)
                    {
                        std::cerr << "WARNING: nom_state.xp_dot too small." << std::endl;
                    }
                }
                else if (real_state.xp_dot <= 1e-1)
                {
                    u << 0, -5;
                    if (VERBOSE)
                    {
                        std::cerr << "WARNING: real_state.xp_dot too small." << std::endl;
                    }
                }
                else
                {
                    // constants
                    //double a = 1.41, b = 1.576, mu = 0.5, Fzf = 21940.0/2, Fzr = 21940.0/2;
                    //double cf = 65000.0, cr = 65000.0, m = 2194.0, Iz = 4770.0;

		    double a = 1.68, b = 1.715, mu = 3.4812e+05, Fzf = 21940.0/2, Fzr = 21940.0/2;    
		    double cf = 3.4812e+05, cr = 3.5537e+05, m = 9840.0, Iz = 41340.0;
                    double psi_dot_com = 0.0; 
 
                    //feedback states:
                    Eigen::Vector2d L_f_output, L_f_f_output, p, p_dot;
                    Eigen::Matrix2d L_g_f_output;
                    double xp_dot, yp_dot, psi_dot, epsi, ey, s;
                    xp_dot = real_state.xp_dot;
                    yp_dot = real_state.yp_dot;
                    psi_dot = real_state.psi_dot;
                    epsi = real_state.epsi;
                    ey= real_state.ey;
                    s = real_state.s;
                    L_f_output <<  yp_dot*cos(epsi) + xp_dot*sin(epsi),
                             xp_dot*cos(epsi) - yp_dot*sin(epsi);
                    L_g_f_output << (2*cf*cos(epsi))/m, sin(epsi),
                            -(2*cf*sin(epsi))/m, cos(epsi);

                    L_f_f_output << (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi),
         sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi);
                    p << ey, s;
                    p_dot = L_f_output;

                    //nominal states:
                    Eigen::Vector2d L_f_output_nom, L_f_f_output_nom, p_nom, p_dot_nom;
                    Eigen::Matrix2d L_g_f_output_nom;
                    xp_dot = nom_state.xp_dot;
                    yp_dot = nom_state.yp_dot;
                    psi_dot = nom_state.psi_dot;
                    epsi = nom_state.epsi;
                    ey= nom_state.ey;
                    s = nom_state.s;
                    L_f_output_nom <<  yp_dot*cos(epsi) + xp_dot*sin(epsi),
                                    xp_dot*cos(epsi) - yp_dot*sin(epsi);
                    L_g_f_output_nom << (2*cf*cos(epsi))/m, sin(epsi),
                                       -(2*cf*sin(epsi))/m, cos(epsi);
                    L_f_f_output_nom << (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi),
							sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi);
                    p_nom << ey, s;
                    p_dot_nom = L_f_output_nom;
 
                    Eigen::Vector2d p_err, p_err_dot;
                    p_err << real_state.ey - nom_state.ey, real_state.s - nom_state.s;
                    p_err_dot = L_f_output - L_f_output_nom;

                    Eigen::Matrix2d k1, k2;
                    k1 << 2.0, 0.0, 0.0, 2.0;
                    k2 <<  2 * 1.414* sqrt(k1(0)), 0.0, 0.0,  2 * 1.414 *sqrt(k1(3));
                    u = -k1 * p_err - k2 * p_err_dot + nom_u;
                    //u= L_g_f_output.inverse()*( -k1*p_err -k2*p_err_dot - L_f_f_output + L_g_f_output_nom*nom_u + L_f_f_output_nom);

                    // The following line is in the original .m file as an alternative output formula
                    // u= L_g_f_output\( -k1*p_err -k2*p_err_dot - L_f_f_output + L_g_f_output_nom*u_nom + L_f_f_output_nom);

                    Eigen::VectorXd xi_err(6); 
                    xi_err<< real_state.xp_dot-nom_state.xp_dot, real_state.yp_dot - nom_state.yp_dot,  real_state.psi_dot - nom_state.psi_dot,
                     real_state.epsi - nom_state.epsi, real_state.ey - nom_state.ey,  real_state.s - nom_state.s; 
                    Eigen::MatrixXd K_state(2,6);
                    K_state <<  -0.0000,   1.7535,    2.3066,   26.2826,    5.4772,   -0.0000,
                                4.5776,    0.0000,   -0.0000,   -0.0000,   -0.0000,   5.4772;  //lqr
                    K_state << 0.0000,    1.6340,    3.6865,   26.1752,    5.4772,    0.0000,
                               4.5776,   -0.0000,    0.0000,   -0.0000,   -0.0000,    5.4772;  //lqr
                    Eigen::Matrix2d k_scale;
                    //k_scale << 1, 0, 0, 1; 
                    k_scale << k_scale_steer, 0, 0, k_scale_acc;
                    u = -k_scale*K_state * xi_err + nom_u;  //more stable 
                    //u = -k_scale*K_state * xi_err;            
                                        
                    //20190216: if the input to the actual system is speed: 
                    u(1) = -0.5 * p_err(1)  + nom_state.xp_dot; 
                    u(1) = -0.5 * (real_state.s - nom_state.s)  + nom_state.xp_dot;   
                    //u(1) = 10;  u(0) = 0; //tune

                    //input is the speed and steering angle: 
                    Eigen::VectorXd xi_err_speed(5); 
                    xi_err_speed << real_state.yp_dot - nom_state.yp_dot,  real_state.psi_dot - nom_state.psi_dot,
                     real_state.epsi - nom_state.epsi, real_state.ey - nom_state.ey,  real_state.s - nom_state.s; 
                    Eigen::MatrixXd K_state_speed(2,5);
                    K_state_speed << 3.6009,   -6.0016,  -17.6119,  -1.7321,    0.0000,
                                     0.0000,   -0.0000,   -0.0000,   -0.0000,   1.7321;  //lqr
                    K_state_speed << 10.5182,  -17.8675,  -54.7332,   -5.4772,  -0.0000,
                                     -0.0000,    0.0000,    0.0000,    0.0000,   5.4772; 
                    //K_state_speed << 4.6047,   -0.5203,  -25.2560,   -5.4772,   -0.0000,
                     //  -0.0000,   0.0000,    0.0000,    0.0000,   5.4772; 

                    Eigen::Matrix2d k_scale_speed;
                    k_scale_speed << 1, 0, 0, 1; 
                    nom_u(1) = nom_state.xp_dot; 
                    Eigen::Vector2d u_tmp;
                    u_tmp = -k_scale_speed*K_state_speed * xi_err_speed + nom_u;  //more stable 
                    //u(1) = u_tmp(1);
                }
                /*opendlv::proxy::PedalPositionRequest pprMsg;
                opendlv::proxy::GroundSteeringRequest gsrMsg;
                pprMsg.position((float)u(0));
                gsrMsg.groundSteering((float)u(1));
                od4_2.send(pprMsg);
                od4_2.send(gsrMsg);*/

                // internal::nomU msgActualu;  //the control variable 
                opendlv::proxy::GroundSpeedRequest msgSpeed;
                opendlv::proxy::GroundSteeringRequest msgSteering;

                if(flag_nomu |  flag_nomstate | flag_actstate) {
                    /*msgActualu.steer(0);
		    msgActualu.acc(0);
                    msgActualu.speed(16); //20190216: if the input to the actual system is speed:
                    */
                    msgSpeed.groundSpeed(16);
                    msgSteering.groundSteering(0);
                }
                else{
                    /*msgActualu.steer(u(0));
		    msgActualu.acc(u(1));
                    msgActualu.speed(u(1));  //20190216: if the input to the actual system is speed:
                    */
                    msgSpeed.groundSpeed(u(1));
                    msgSteering.groundSteering(u(0));
                }
                //od4_2.send(msgSpeed);
                //od4_2.send(msgSteering);

                std::ofstream txt2("/tmp/data_msg_actual_u.txt", std::ios::out | std::ios::app);
                if (txt2.is_open())
                {
                    //the clock is not accurate: 
                    txt2 << ((double)clock())/CLOCKS_PER_SEC << '\t' << msgSpeed.groundSpeed() << '\t' << msgSteering.groundSteering() << '\n';
                    txt2.close();
                }
                else std::cerr << "WARNING: Unable to save data into the file <data_msg_nom_u.txt>." << std::endl;

                if (VERBOSE)
                {
                    std::cout << "Request messages sent:" << std::endl << "[" << u(0) << ", " << u(1) << "]" << std::endl;
                }
                return false;
            } // end of inner lambda function
        }; // end of sendMsg
        od4_2.timeTrigger(FREQ,sendMsg);
    }
    return 0;
}

