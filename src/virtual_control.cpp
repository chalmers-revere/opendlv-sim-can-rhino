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
#include <fstream>
#include <thread>
#include <ctime>

#include "cluon-complete.hpp"
#include "opendlv-standard-message-set.hpp"

#include "data_structure.hpp"

int32_t main(int32_t argc, char *argv[])
{
    int  flag_ini = 0;    //initialization flag, 
    auto commandlineArguments = cluon::getCommandlineArguments(argc, argv);

    uint16_t CID = 113;
    if (0 == commandlineArguments.count("cid"))
    {
        std::cerr << "WARNING: No cid assigned, using 113 by default ('--cid=113')." << std::endl;
    }
    else CID = std::stoi(commandlineArguments["cid"]);

    uint32_t FREQ = 50;
    if (0 == commandlineArguments.count("freq"))
    {
        std::cerr << "WARNING: No frequency assigned, using 50 by default ('--freq=50')." << std::endl;
    }
    else FREQ = std::stoi(commandlineArguments["freq"]);

    double v_ref = 16.0;
    if (0 == commandlineArguments.count("v_ref"))
    {
        std::cerr << "WARNING: No v_ref assigned, using 16.0 by default ('--v_ref=16')." << std::endl;
    }
    else v_ref = std::stoi(commandlineArguments["v_ref"]);
 
    int no_obs = 3;
    if (0 == commandlineArguments.count("no_obs"))
    {
        std::cerr << "WARNING: No no_obs assigned, using 16.0 by default ('--no_obs=3')." << std::endl;
    }
    else no_obs = std::stoi(commandlineArguments["no_obs"]);

    double pos_ycenter_ob = -1.2;
    if (0 == commandlineArguments.count("pos_ycenter_ob"))
    {
        std::cerr << "WARNING: No pos_ycenter_ob assigned, using 16.0 by default ('--pos_ycenter_ob=-1.2')." << std::endl;
    }
    else pos_ycenter_ob = std::stoi(commandlineArguments["pos_ycenter_ob"]);

    bool const VERBOSE{commandlineArguments.count("verbose") != 0};

    FB_state nom_state(v_ref, 0, 0, 0, 0, 0, 0, 0);
    if (VERBOSE) 
        std::cout << "Nom_state initialised." << std::endl;

    auto currentTime = std::to_string(cluon::time::toMicroseconds(cluon::time::now()) / 1000 / 60 ); // resolution to minutes
    auto filename = "/tmp/data_traj_ob_" + currentTime + ".txt";
    auto filename2 = "/tmp/data_msg_nom_u_" + currentTime + ".txt";
    auto filename3 = "/tmp/data_ref_traj_" + currentTime + ".txt";

    Global_variables gl(v_ref);
    gl.no_ob = no_obs;
    gl.isVerbose = VERBOSE;
    gl.generate_init_ob(pos_ycenter_ob);
    gl.currentTime = currentTime;

    // Data saving into txt file
    // number of rows = number of obstacles
    // each row contains the following data, seperated by tab
    // pos_x  pos_y  vel_x  vel_y  acc_x  acc_y  radius

   /* std::ofstream txt("/tmp/data_traj_ob.txt", std::ios::out);
    if (txt.is_open())
    {
        for (uint8_t i = 0; i < gl.no_ob; ++i)
        {
            Obstacle curr = gl.traj_ob[i];
            txt << curr.pos_x << '\t' << curr.pos_y << '\t' << curr.vel_x << '\t' << curr.vel_y << '\t' 
                << curr.acc_x << '\t' << curr.acc_y << '\t' << curr.radius << '\n';
        }
        txt.close();
    }
    else std::cerr << "WARNING: Unable to save data into the file <data_traj_ob.txt>." << std::endl;  */


    if (VERBOSE)
    {
        std::cout << "Obstacles generated." << std::endl;
    }

    cluon::OD4Session od4(CID, [&nom_state, &VERBOSE, &flag_ini, &gl](cluon::data::Envelope &&env) noexcept {
        // std::cout << "Received message:" << env.dataType() << std::endl;

        if (env.dataType() == internal::nomState::ID())
        {
            internal::nomState received = cluon::extractMessage<internal::nomState>(std::move(env));
            // nom_state = FB_state(received.xp_dot(), received.yp_dot(), received.psi_dot(), 
            //    received.epsi(), received.ey(), received.s(), received.steer(), received.acc());
            nom_state.xp_dot = received.xp_dot();
            nom_state.yp_dot = received.yp_dot();
            nom_state.psi_dot = received.psi_dot();
            nom_state.epsi = received.epsi();
            nom_state.ey = received.ey();
            nom_state.s = received.s();
            nom_state.steer = received.steer();
            nom_state.acc = received.acc();
            //std::cout << "testing flag_ini " <<   flag_ini  << std::endl;
            if(flag_ini == 0) {
            //initialization the reference trajectory: 
               gl.trajd[2] << 0, 0, 0;  //tra_com_ddot
               gl.trajd[1] << 0, 0, gl.v_ref;  //tra_com_dot
               gl.trajd[0] << 0, 0, nom_state.s;   //tra_com
            }

            if (VERBOSE)
            {
               std::cout << "New nom_state received:" << std::endl;
               nom_state.print();  
               std::cout << "flag_ini: " << flag_ini << std::endl;
            }
            flag_ini++; if (flag_ini == 127) flag_ini = 10;
        }
    });

    if (0 == od4.isRunning())
    {
        std::cerr << "ERROR: No OD4 running!!!" << std::endl;
        return -1;
    }
    while (od4.isRunning()) //should run after initialization 
    {
        auto sendMsg{[&od4, &nom_state, &gl, &VERBOSE, &FREQ, &flag_ini, &filename, &filename2, &filename3]() -> bool
            {
                // update position of obstacles
                gl.ob_traj(false); 

                // update trajd
                gl.traj_gen(nom_state, FREQ);

                //20190108:
		//gl.trajd[0](2) = gl.trajd[0](2) + 16.0/FREQ;  

                // run the solver
                Output_safety correct = safety_certificate_complex(nom_state, gl);

		//tunning, why always the same? 20190103
		if (VERBOSE) {
                  std::cout << "current state in virtual control:" << std::endl;
                  nom_state.print(); 
                }

                gl.nosolution = !(correct.hasSolution);

                internal::nomU msgNomU;

                if (flag_ini >= 1){
                //  notice the order:                 
                msgNomU.steer(correct.x(0));
		msgNomU.acc(correct.x(1));
                }
                else{
                msgNomU.steer(0);
		msgNomU.acc(0);
                }

               // msgNomU.steer(0);
	       // msgNomU.acc(1);
                od4.send(msgNomU);

		
                if (VERBOSE)
                {	            
                    std::cout << "Message nomU sent: " << std::endl << "[" << msgNomU.acc() << ", " << msgNomU.steer() << "]" << std::endl;
                }

                // Data saving into txt file
                // each row contains the following data, seperated by tab
                // time msgNomU.acc  msgNomU.steer
                std::ofstream txt2(filename2, std::ios::out | std::ios::app);
                if (txt2.is_open())
                {
                    txt2 << ((double)clock())/CLOCKS_PER_SEC << '\t' << msgNomU.acc() << '\t' << msgNomU.steer() << '\n';
                    txt2.close();
                }
                else std::cerr << "WARNING: Unable to save data into the file <data_msg_nom_u.txt>." << std::endl;

                std::ofstream txt3(filename3, std::ios::out | std::ios::app);
                if (txt3.is_open())
                {
                    txt3 << gl.trajd[0](0) << '\t' << gl.trajd[0](1)  << '\t' << gl.trajd[0](2)   << '\t' << gl.trajd[1](0) << '\t' << gl.trajd[1](1)  << '\t' << gl.trajd[1](2) << '\n';
                    txt3.close();
                }
                else std::cerr << "WARNING: Unable to save data into the file <data_ref_traj.txt>." << std::endl;

                std::ofstream txt(filename, std::ios::out | std::ios::app );
                if (txt.is_open())
                {
                for (uint8_t i = 0; i < gl.no_ob; ++i){
                    Obstacle curr = gl.traj_ob[i];
                    txt << curr.pos_x << '\t';                    
                }
                txt <<  '\n';
                for (uint8_t i = 0; i < gl.no_ob; ++i){
                    Obstacle curr = gl.traj_ob[i];
                    txt << curr.pos_y << '\t';                    
                }
                txt <<  '\n';
                for (uint8_t i = 0; i < gl.no_ob; ++i){
                    Obstacle curr = gl.traj_ob[i];
                    txt << curr.radius << '\t';                    
                }
                txt <<  '\n';
                txt.close();
                }
                else std::cerr << "WARNING: Unable to save data into the file <data_traj_ob.txt>." << std::endl;
                return false;
            }
        };
        od4.timeTrigger(FREQ,sendMsg);
    }
    return 0;
}
