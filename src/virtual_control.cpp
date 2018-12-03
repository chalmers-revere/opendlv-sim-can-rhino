#include <cstdint>
#include <chrono>
#include <iostream>
#include <sstream>
#include <thread>

#include "cluon-complete.hpp"

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

    FB_state nom_state(0, 0, 0, 0, 0, 0, 0, 0);
    if (VERBOSE) 
        std::cout << "Nom_state initialised." << std::endl;

    Global_variables gl;

    cluon::OD4Session od4(CID, [](cluon::data::Envelope &&env) noexcept {
        if (env.dataType() == internal::nomState::ID())
        {
            internal::nomState received = cluon::extractMessage<internal::nomState>(std::move(env));
            nom_state = FB_state(received.xp_dot, received.yp_dot, received.epsi, received.ey, received.s, received.steer, received.acc);
            if (VERBOSE)
            {
                std::cout << "New nom_state received:" << std::endl;
                nom_state.print(); 
            }
        }

        if (env.dataType() == internal::nomState::ID())
    });

    if (0 == od4.isRunning())
    {
        std::cerr << "ERROR: No OD4 running!!!" << std::endl;
        return -1;
    }
    else
    {
        auto sendMsg{[&od4, &nom_state, &gl]() -> bool
            {
                //gl.scale unused
                Output_safety correct = safety_certificate_complex(nom_state, gl);
                gl.nosolution = !(correct.hasSolution);

                internal::nomU msgNomU;
                if (gl.nosolution)
                {
                    std::cerr << "WARNING: No solution detected from solver!!" << std::endl;
                    msgNomU.acc = 0.0;
                    msgNomU.steer = 0.0;
                }
                else
                {
                    msgNomU.acc = correct.x(0);
                    msgNomU.steer = correct.x(1);
                }
                od4.send(msgNomU);
                if (VERBOSE)
                {
                    std::cout << "Message nomU sent: [" << msgNomU.acc << ", " << msgNomU.steer << "]" << std::endl;
                }
                return false;
            }
        };
        od4.timeTrigger(FREQ,sendMsg);
    }

    return 0;
}

/* Eigen::Vector2d virtual_control(FB_state state, Global_variables& gl)
{
    Eigen::Vector2d u;
    if (0 == gl.scale)
    {
        // "horizon = 1" unused
        Output_safety correct = safety_certificate_complex(state, gl); 
        gl.nosolution = !(correct.hasSolution);
        u = correct.x;
        gl.u_global = u;
    }
    else u = gl.u_global;

    gl.scale = (gl.scale < 20) ? gl.scale + 1 : 0;

    return u;
} */