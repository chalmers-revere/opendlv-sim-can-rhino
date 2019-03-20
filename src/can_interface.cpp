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

// GPS converter
#include "WGS84toCartesian.hpp"

#include <cstdint>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <thread>
#include "dynamics_vehicle.h"
#include <math.h>

#include <cmath>
#include <Eigen/Dense>

int32_t main(int32_t argc, char **argv) 
{ 
    bool flag_ini_position{false};
    bool flag_if_dataTriggered{false};
   
    // int32_t retCode{0};

    auto commandlineArguments = cluon::getCommandlineArguments(argc, argv);
    if ((0 == commandlineArguments.count("cid")) 
        || (0 == commandlineArguments.count("gps-id")) 
        || (0 == commandlineArguments.count("can-id")) 
        || (0 == (commandlineArguments.count("freq") ^ commandlineArguments.count("dt")))) 
    {
        std::cerr << argv[0] << "CAN interface." << std::endl;
        std::cerr << "Usage:   " << argv[0] << " --cid=<OpenDaVINCI session> --gps-id=<Sender ID of GPS message> --can-id=<Sender ID of CAN message> [--freq=<Frequency>] [--dt<if sending in timeTriggered way>] [--verbose]" << std::endl;
        /* TODO: Verify if necessary to separate inout and output ID */
        std::cerr << "Example: " << argv[0] << " --cid=114 --gps-id=0 --can-id=1 [--freq=50 or --dt] [--verbose]" << std::endl;
        // std::cerr << "(The second OD4Session is for vehicle state overriding from external source.)" << std::endl;
        return -1;
    } else {
        const bool VERBOSE{commandlineArguments.count("verbose") != 0};
        uint32_t FREQ;
        if (commandlineArguments.count("dt") != 0)
        {
            flag_if_dataTriggered = true;
            FREQ = 0;
        }
        else
        {
            FREQ = static_cast<uint32_t>(std::stoi(commandlineArguments["freq"]));
        }

        uint16_t argCID = static_cast<uint16_t>(std::stoi(commandlineArguments["cid"]));
        uint16_t argGPS_ID = static_cast<uint16_t>(std::stoi(commandlineArguments["gps-id"]));
        uint16_t argCAN_ID = static_cast<uint16_t>(std::stoi(commandlineArguments["can-id"]));

        std::shared_ptr<cluon::OD4Session> od4_2 = std::shared_ptr<cluon::OD4Session>(new cluon::OD4Session(argCID));
        /*"od4_2" means "external channel" by default*/

        typedef struct {
            double init_latitude{0.0};
            double init_longitude{0.0};
            double init_heading{0.0};
            double s{0.0};
            double ey{0.0};
            double omega_body[3]{};
            double v_body{16.0};
            float heading{0.0};
            Eigen::Vector3d v_ned; 
            Eigen::Vector3d v_body_3d; 
            Eigen::Vector3d Eulerangle; 
        }All_vars;

        All_vars var;

        while (od4_2->isRunning())
        {

            auto Reading_GPS{[&var, &VERBOSE, &flag_ini_position, &argGPS_ID](cluon::data::Envelope &&env)
                {
                    if (env.senderStamp() != argGPS_ID)
                    {
                        return;
                    }
                    opendlv::logic::sensation::Geolocation reading_wgs84 = cluon::extractMessage<opendlv::logic::sensation::Geolocation>(std::move(env));
                    if (!flag_ini_position) // first reading of GPS i.e. no ref point set yet
                    {
                        var.init_latitude = reading_wgs84.latitude();
                        var.init_longitude = reading_wgs84.longitude();
                        var.init_heading = reading_wgs84.heading();
                        flag_ini_position = true;
                        if (VERBOSE) 
                            std::cout << "Init position ref_point set." << std::endl;
                    }
                    else
                    {
                        std::array<double, 2> ref_point{var.init_latitude, var.init_longitude};
                        std::array<double, 2> pos_wgs84{reading_wgs84.latitude(), reading_wgs84.longitude()};
                        std::array<double, 2> pos_cartesian{wgs84::toCartesian(ref_point, pos_wgs84)};

                        var.s = pos_cartesian[0];
                        var.ey = pos_cartesian[1];
                        var.heading = reading_wgs84.heading();
                        if (VERBOSE)
                        {
                            std::cout << "Received GPS reading, WGS84: [" << pos_wgs84[0] << ", " << pos_wgs84[1] << "], ";
                            std::cout << "Cartesian: [" << pos_cartesian[0] << ", " << pos_cartesian[1] << "]." << std::endl;
                            std::cout << "Heading: " << reading_wgs84.heading() << std::endl;
                        }
                    }
                }
            }; // end of Reading_GPS

            auto Reading_AngularVelocity{[&var, &VERBOSE, &argGPS_ID](cluon::data::Envelope &&env)
                {
                    if (env.senderStamp() != argGPS_ID)
                    {
                        return;
                    }
                    opendlv::proxy::AngularVelocityReading reading_ang_vel = cluon::extractMessage<opendlv::proxy::AngularVelocityReading>(std::move(env));
                    var.omega_body[0] = reading_ang_vel.angularVelocityX();
                    var.omega_body[1] = reading_ang_vel.angularVelocityY();
                    var.omega_body[2] = reading_ang_vel.angularVelocityZ();
                    if (VERBOSE)
                    {
                        std::cout << "Received angular velocity reading: [" 
                            << reading_ang_vel.angularVelocityX() << ", " 
                            << reading_ang_vel.angularVelocityY() << ", " 
                            << reading_ang_vel.angularVelocityZ() << "]." << std::endl;
                    }
                }
            }; // end of Reading_AngularVelocity

            auto Reading_YawRate{[&var, &VERBOSE, &argCAN_ID](cluon::data::Envelope &&env)
                {
                    if (env.senderStamp() != argCAN_ID)
                    {
                        return;
                    }
                    opendlv::proxy::rhino::VehicleState reading_yr = cluon::extractMessage<opendlv::proxy::rhino::VehicleState>(std::move(env));

                    // var.omega_body[0] = reading_ang_vel.angularVelocityX();
                    // var.omega_body[1] = reading_ang_vel.angularVelocityY();
                    var.omega_body[2] = reading_yr.yawRate();
                    if (VERBOSE)
                    {
                        std::cout << "Received yawRate reading (from CAN): " 
                            << reading_yr.yawRate() << std::endl;
                    }
                }
            }; // end of Reading_YawRate

            auto Reading_GroundSpeed{[&var, &VERBOSE, &argCAN_ID](cluon::data::Envelope &&env)
                {
                    if (env.senderStamp() != argCAN_ID)
                    {
                        return;
                    }
                    opendlv::proxy::GroundSpeedReading reading_gs = cluon::extractMessage<opendlv::proxy::GroundSpeedReading>(std::move(env));
                    // TODO: verify the following line (regarding frames and potential need of conversion
                    var.v_body = reading_gs.groundSpeed();

                    double phi= 0;
                    double theta= 0;
                    double psi= -var.heading;  //down: z 
                //rotation matrix of body frame relating to inertial frame:
                    Eigen::Matrix3d R_rot; 
                    R_rot << cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), 
                           cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), 
                          -sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta); 
                    //var.v_ned; 
                    //Eigen::Vector3d v_body;  
                    var.v_body_3d = R_rot.adjoint() * var.v_ned;  //front x, right y, down z
                    var.v_body_3d(1) =  -var.v_body_3d(1);   //left y
                    var.v_body_3d(2) =  -var.v_body_3d(2);   //up z

                    if (VERBOSE)
                    {
                        std::cout << "Received ground speed reading: " << reading_gs.groundSpeed() << std::endl;
                    }
                }
            }; // end of Reading_GroundSpeed

            if (od4_2->isRunning())
            {
                od4_2->dataTrigger(opendlv::logic::sensation::Geolocation::ID(), Reading_GPS);
                od4_2->dataTrigger(opendlv::proxy::AngularVelocityReading::ID(), Reading_AngularVelocity);
                od4_2->dataTrigger(opendlv::proxy::GroundSpeedReading::ID(), Reading_GroundSpeed);
                od4_2->dataTrigger(opendlv::proxy::rhino::VehicleState::ID(), Reading_YawRate);
            }

            auto Output{[&var, &VERBOSE, od4_2, &argCAN_ID]() -> bool
                {
                    opendlv::sim::Frame msg;
                    double psi_ini = var.init_heading + 3.14159265/2; 
                    double s = cos(psi_ini) * var.s + sin(psi_ini) * var.ey; 
                    double ey = -sin(psi_ini) * var.s + cos(psi_ini) * var.ey;
                    double psi =  var.heading - var.init_heading; 
                    msg.x((float)s);
                    msg.y((float)ey);
                    msg.yaw((float)psi);
                    msg.z(0.0);
                    msg.roll(0.0);
                    msg.pitch(0.0);
                    od4_2->send(msg, cluon::time::now(), argCAN_ID);

                    opendlv::sim::KinematicState msg_2;
                    msg_2.vx((float)var.v_body_3d(0));
                    msg_2.vy((float)var.v_body_3d(1));
                    msg_2.yawRate((float)var.omega_body[2]);
                    msg_2.vz(0.0);
                    msg_2.rollRate(0.0);
                    msg_2.pitchRate(0.0);
                    od4_2->send(msg_2, cluon::time::now(), argCAN_ID);

                    internal::nomState msg_3;
                    msg_3.xp_dot((float)var.v_body_3d(0));
                    msg_3.yp_dot((float)var.v_body_3d(1));
                    msg_3.psi_dot(var.omega_body[2]);
                    msg_3.epsi(psi);
                    msg_3.ey(ey);
                    msg_3.s(s);
                    msg_3.steer(0);
                    msg_3.acc(0); 
                    od4_2->send(msg_3, cluon::time::now(), argCAN_ID);

                    if (VERBOSE)
                    {
                        std::cout << "Position and Kinematic state messages sent." << std::endl;
                    }

                        // Data saving into txt file (if "save_file" indicated)
                        // number of rows = length of time 
                        // each row contains the following data, seperated by tab: 
                        // time nomStateMsg(all attributes)
                    {
                        std::ofstream txt("/tmp/data_model_state.txt", std::ios::out | std::ios::app);
                        if (txt.is_open())
                        {
                                txt << 0 << '\t'
                                << msg_3.xp_dot() << '\t' << msg_3.yp_dot() << '\t' 
                                << msg_3.psi_dot() << '\t' << msg_3.epsi() << '\t' 
                                << msg_3.ey() << '\t' << msg_3.s() << '\t' 
                               << msg_3.steer() << '\t' << msg_3.acc() << '\n';
                                txt.close();
                        }
                        else std::cerr << "WARNING: Unable to save data into the file <" << "/tmp/data_model_state.txt" << ">." << std::endl;
                     }

                    return false;
                }
            };

            auto Output_dataTriggered{[&var, &VERBOSE, od4_2, &argCAN_ID](cluon::data::Envelope &&env) -> bool
                {
                    internal::nomU msgU = cluon::extractMessage<internal::nomU>(std::move(env));
                    if (VERBOSE)
                    {
                        std::cout << "Received nomU, sending position and kinematic state messages ...";
                    }
                    
                    opendlv::sim::Frame msg;
                    double psi_ini = var.init_heading + 3.14159265/2; 
                    double s = cos(psi_ini) * var.s + sin(psi_ini) * var.ey; 
                    double ey = -sin(psi_ini) * var.s + cos(psi_ini) * var.ey;
                    double psi =  var.heading - var.init_heading; 
                    msg.x((float)s);
                    msg.y((float)ey);
                    msg.yaw((float)psi);
                    msg.z(0.0);
                    msg.roll(0.0);
                    msg.pitch(0.0);
                    //od4_2->send(msg, cluon::time::now(), argCAN_ID);

                    opendlv::sim::KinematicState msg_2;
                    msg_2.vx((float)var.v_body_3d(0));
                    msg_2.vy((float)var.v_body_3d(1));
                    msg_2.yawRate((float)var.omega_body[2]);
                    msg_2.vz(0.0);
                    msg_2.rollRate(0.0);
                    msg_2.pitchRate(0.0);
                    //od4_2->send(msg_2, cluon::time::now(), argCAN_ID);

                    internal::nomState msg_3;
                    msg_3.xp_dot((float)var.v_body_3d(0));
                    msg_3.yp_dot((float)var.v_body_3d(1));
                    msg_3.psi_dot(var.omega_body[2]);
                    msg_3.epsi(psi);
                    msg_3.ey(ey);
                    msg_3.s(s);
                    msg_3.steer(0);
                    msg_3.acc(0); 
                    od4_2->send(msg_3, cluon::time::now(), argCAN_ID);

                    if (VERBOSE)
                    {
                        std::cout << "Sent." << std::endl;
                    }

                    return false;
                }
            };

            if (od4_2->isRunning())
            {
                if (flag_if_dataTriggered)
                {
                    od4_2->dataTrigger(internal::nomU::ID(), Output_dataTriggered);
                }
                else
                {
                    od4_2->timeTrigger(FREQ, Output);
                }
            }

        } // end of while
    } // end of else
    return 1;
}
