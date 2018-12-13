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
    else CID = std::stoi(commandlineArguments["freq"]);

    bool const VERBOSE{commandlineArguments.count("verbose") != 0};

    FB_state nom_state(4.0, 0, 0, 0, 0, 0, 0, 0);
    if (VERBOSE) 
        std::cout << "Nom_state initialised." << std::endl;

    Global_variables gl;
    gl.isVerbose = VERBOSE;

    gl.generate_init_ob();

    if (VERBOSE)
    {
        std::cout << "Obstacles generated." << std::endl;
    }

    cluon::OD4Session od4(CID, [&nom_state, &VERBOSE](cluon::data::Envelope &&env) noexcept {
        if (env.dataType() == internal::nomState::ID())
        {
            internal::nomState received = cluon::extractMessage<internal::nomState>(std::move(env));
            nom_state = FB_state(received.xp_dot(), received.yp_dot(), received.psi_dot(), received.epsi(), received.ey(), received.s(), received.steer(), received.acc());
            if (VERBOSE)
            {
                std::cout << "New nom_state received:" << std::endl;
                nom_state.print(); 
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
        auto sendMsg{[&od4, &nom_state, &gl, &VERBOSE]() -> bool
            {
                // update position of obstacles
                gl.ob_traj(true); 
                /*if (VERBOSE)
                {
                    std::cout << "Obstacle status:" << std::endl;
                    for (uint16_t i = 0; i < gl.no_ob; ++i)
                    {
                        std::cout << "No. " << i << " ";
                        gl.traj_ob[i].print();
                    }
                }*/

                // update trajd
                gl.traj_gen(nom_state);

                // run the solver
                Output_safety correct = safety_certificate_complex(nom_state, gl);
                gl.nosolution = !(correct.hasSolution);

                internal::nomU msgNomU;
                if (gl.nosolution)
                {
                    std::cerr << "WARNING: No solution detected from solver!!" << std::endl;
                    msgNomU.acc(0.0);
                    msgNomU.steer(0.0);
                }
                else
                {
                    msgNomU.acc(correct.x(0));
                    msgNomU.steer(correct.x(1));
                }
                od4.send(msgNomU);
                if (VERBOSE)
                {
                    std::cout << "Message nomU sent: " << std::endl << "[" << msgNomU.acc() << ", " << msgNomU.steer() << "]" << std::endl;
                }
                return false;
            }
        };
        od4.timeTrigger(FREQ,sendMsg);
    }
    return 0;
}