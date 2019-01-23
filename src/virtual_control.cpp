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

    bool const VERBOSE{commandlineArguments.count("verbose") != 0};

    FB_state nom_state(6.0, 0, 0, 0, 0, 0, 0, 0);
    if (VERBOSE) 
        std::cout << "Nom_state initialised." << std::endl;

    Global_variables gl;
    gl.isVerbose = VERBOSE;

    gl.generate_init_ob();

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

    cluon::OD4Session od4(CID, [&nom_state, &VERBOSE](cluon::data::Envelope &&env) noexcept {
        std::cout << "Received message:" << env.dataType() << std::endl;

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

            if (VERBOSE)
            {

               //20190103, where is the nom_state? cannot see it. 
               // std::cout << "New nom_state received:" << std::endl;
               // nom_state.print(); 
            }
        }
    });

    if (0 == od4.isRunning())
    {
        std::cerr << "ERROR: No OD4 running!!!" << std::endl;
        return -1;
    }
    while (od4.isRunning())
    {
        auto sendMsg{[&od4, &nom_state, &gl, &VERBOSE, &FREQ]() -> bool
            {
                // update position of obstacles
                gl.ob_traj(false); 

                // update trajd
                //gl.traj_gen(nom_state);

                //20190108:
		gl.trajd[0](2) = gl.trajd[0](2) + 12.0/FREQ;  

                // run the solver
                Output_safety correct = safety_certificate_complex(nom_state, gl);

		//tunning, why always the same? 20190103
		if (VERBOSE) {
                  std::cout << "current state in virtual control:" << std::endl;
                  nom_state.print(); 
                }

                gl.nosolution = !(correct.hasSolution);

                internal::nomU msgNomU;
                //  notice the order:                 
                msgNomU.steer(correct.x(0));
		msgNomU.acc(correct.x(1));

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
                std::ofstream txt2("/tmp/data_msg_nom_u.txt", std::ios::out | std::ios::app);
                if (txt2.is_open())
                {
                    txt2 << ((double)clock())/CLOCKS_PER_SEC << '\t' << msgNomU.acc() << '\t' << msgNomU.steer() << '\n';
                    txt2.close();
                }
                else std::cerr << "WARNING: Unable to save data into the file <data_msg_nom_u.txt>." << std::endl;

                std::ofstream txt("/tmp/data_traj_ob.txt", std::ios::out | std::ios::app );
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
