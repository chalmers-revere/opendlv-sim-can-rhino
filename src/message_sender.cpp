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
#include <fstream>
#include <sstream>
#include <thread>

#include "cluon-complete.hpp"
#include "opendlv-standard-message-set.hpp"

// #include "data_structure.hpp"

int32_t main(int32_t argc, char *argv[])
{
    auto commandlineArguments = cluon::getCommandlineArguments(argc, argv);

    const bool ifRead{commandlineArguments.count("read_file") != 0};
    if (!ifRead)
    {
        std::cerr << "ERROR: No data file indicated." << std::endl;
        std::cerr << "\tExample: --read_file=/tmp/data_filename" << std::endl;
        return -1;
    }

    auto filename{commandlineArguments["read_file"]};

    std::ifstream txt(filename, std::ios::in);
    if ( !(txt.is_open()) )
    {
        std::cerr << "ERROR: Data file not available." << std::endl;
        return -2;
    }

    if (0 == commandlineArguments.count("msg_id"))
    {
        std::cerr << "ERROR: No message type indicated." << std::endl;
        std::cerr << "\tExample: --msg_id=3001\n\t(3001 for internal.nomU and 3002 for internal.nomState)" << std::endl;
        return -3;
    }
    uint16_t msgID = std::stoi(commandlineArguments["msg_id"]);
    std::cout << "NOTE: Structure of data in the txt file is NOT checked automatically for the time being." << std::endl;
    std::cout << "Please ensure the alignment of the data in the file." << std::endl;

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

    std::shared_ptr<cluon::OD4Session> od4 = std::shared_ptr<cluon::OD4Session>(new cluon::OD4Session(CID));

    while (od4->isRunning() && !(txt.fail()))
    {
        auto sendMsg{[&od4, &VERBOSE, &txt, &msgID]() -> bool
            {
                if (msgID == internal::nomU::ID())
                {
                    internal::nomU msgNomU;
                    double time; // unused for broadcast
                    double acc_read, steer_read;
                    txt >> time >> acc_read >> steer_read;
                    msgNomU.acc(acc_read);
                    msgNomU.steer(steer_read);
                    od4->send(msgNomU);
                    if (VERBOSE)
                    {
                        std::cout << "Message nomU sent: " << std::endl << "[" << msgNomU.acc() << ", " << msgNomU.steer() << "]" << std::endl;
                    }
                }
                else if (msgID == internal::nomState::ID())
                {
                    internal::nomState nomStateMsg;
                    double time; // unused for broadcast
                    double xp_dot_read, yp_dot_read, psi_dot_read, epsi_read, ey_read, s_read, steer_read, acc_read;
                    txt >> time >> xp_dot_read >> yp_dot_read >> psi_dot_read >> epsi_read
                        >> ey_read >> s_read >> steer_read >> acc_read;
                    nomStateMsg.xp_dot(xp_dot_read);
                    nomStateMsg.yp_dot(yp_dot_read);
                    nomStateMsg.psi_dot(psi_dot_read);
                    nomStateMsg.epsi(epsi_read);
                    nomStateMsg.ey(ey_read);
                    nomStateMsg.s(s_read);
                    nomStateMsg.steer(steer_read);
                    nomStateMsg.acc(acc_read);
                    od4->send(nomStateMsg);
                    if (VERBOSE)
                    {
                        std::cout << "Message nomStateMsg sent: " << std::endl;
                        std::cout << "[" << nomStateMsg.xp_dot() << ", " << nomStateMsg.yp_dot() << ", "
                            << nomStateMsg.psi_dot() << ", " << nomStateMsg.epsi() << ", "
                            << nomStateMsg.ey() << ", " << nomStateMsg.s() << ", "
                            << nomStateMsg.steer() << ", " << nomStateMsg.acc() << "]" << std::endl;
                    }
                }
                return false;
            } // end of lambda function
        }; // end of sendMsg definition
        od4->timeTrigger(FREQ,sendMsg);
    }
    return 0;
}

