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



    FB_state nom_state(3.0, 0, 0, 0, 0, 0, 0, 0);
    if (VERBOSE) 
        std::cout << "nom_state initialised." << std::endl;

    FB_state real_state(3.0, 0, 0, 0, 0, 0, 0, 0);

    Eigen::Vector2d nom_u; // Order: [acc, steer]'
    nom_u << 0.0, 0.0;
    if (VERBOSE) 
        std::cout << "nom_u initialised." << std::endl;

    // Global_variables gl;

    cluon::OD4Session od4(CID, [&nom_state, &real_state, &nom_u, &VERBOSE](cluon::data::Envelope &&env) noexcept {
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
        else if (env.dataType() == internal::nomU::ID())
        {
            internal::nomU received = cluon::extractMessage<internal::nomU>(std::move(env));
            nom_u << received.acc(), received.steer();
            if (VERBOSE)
            {
                std::cout << "New nom_u received:" << std::endl;
                std::cout << "[" << nom_u(0) << ", " << nom_u(1) << "]" << std::endl;
            }
        }
        else if (env.dataType() == opendlv::sim::Frame::ID())
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
    });

    if (0 == od4.isRunning())
    {
        std::cerr << "ERROR: No OD4 running!!!" << std::endl;
        return -1;
    }
    while (od4.isRunning())
    {
        auto sendMsg{[&od4, &nom_state, &real_state, &nom_u, &VERBOSE]() -> bool
            {
                Eigen::Vector2d u;
                if (nom_state.xp_dot <= 1e-1)
                {
                    u = nom_u;
                    if (VERBOSE)
                    {
                        std::cerr << "WARNING: nom_state.xp_dot too small." << std::endl;
                    }
                }
                else if (real_state.xp_dot <= 1e-1)
                {
                    u = nom_u;
                    if (VERBOSE)
                    {
                        std::cerr << "WARNING: real_state.xp_dot too small." << std::endl;
                    }
                }
                else
                {
                    // constants
                    double a = 1.41, b = 1.576, mu = 0.5, Fzf = 21940.0/2, Fzr = 21940.0/2;
                    double cf = 65000.0, cr = 65000.0, m = 2194.0, Iz = 4770.0;
                    double psi_dot_com = 0.0, p = Iz / (m * b);

                    Eigen::Vector2d L_f_output, L_f_output_nom;
                    L_f_output << real_state.yp_dot * cos(real_state.epsi) + real_state.xp_dot * sin(real_state.epsi), 
                        real_state.xp_dot * cos(real_state.epsi) - real_state.yp_dot * sin(real_state.epsi); 
                    L_f_output_nom << nom_state.yp_dot * cos(nom_state.epsi) + nom_state.xp_dot * sin(nom_state.epsi), 
                        nom_state.xp_dot * cos(nom_state.epsi) - nom_state.yp_dot * sin(nom_state.epsi); 

                    /*
                    // the following code snippet is in the original .m file, using vehicle states (real_state)
                    L_f_f_output = [ (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi)
                     sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi)];
                    L_g_f_output = [[(2*cf*cos(epsi))/m, sin(epsi)]
                    [ -(2*cf*sin(epsi))/m, cos(epsi)]];
                    p = [ y(5); y(6)];  %output
                    p_dot = L_f_output;  %derivative of output
                    */    
                    /*
                    // the following code snippet is in the original .m file, using nominal states (nom_state)
                    L_f_f_output_nom  = [ (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi)
                     sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi)];
                    L_g_f_output_nom  = [[(2*cf*cos(epsi))/m, sin(epsi)]
                    [ -(2*cf*sin(epsi))/m, cos(epsi)]];
                    p_nom = [y_nom(5); y_nom(6)];  % output, nominal 
                    p_nom_dot = L_f_output_nom;  %derivative of output, nominal 
                    */

                    Eigen::Vector2d p_err, p_err_dot;
                    p_err << real_state.ey - nom_state.ey, real_state.s - nom_state.s;
                    p_err_dot = L_f_output - L_f_output_nom;

                    Eigen::Matrix2d k1, k2;
                    k1 << -1.0, 0.0, 0.0, -3.0;
                    k2 << -2 * 1.414, 0.0, 0.0, -2 * 1.732 * 1.414;
                    u = k1 * p_err + k2 * p_err_dot + nom_u;

                    // The following line is in the original .m file as an alternative output formula
                    // u= L_g_f_output\( -k1*p_err -k2*p_err_dot - L_f_f_output + L_g_f_output_nom*u_nom + L_f_f_output_nom);
                }
                opendlv::proxy::PedalPositionRequest pprMsg;
                opendlv::proxy::GroundSteeringRequest gsrMsg;
                pprMsg.position((float)u(0));
                gsrMsg.groundSteering((float)u(1));
                od4.send(pprMsg);
                od4.send(gsrMsg);
                if (VERBOSE)
                {
                    std::cout << "Request messages sent:" << std::endl << "[" << u(0) << ", " << u(1) << "]" << std::endl;
                }
                return false;
            } // end of inner lambda function
        }; // end of sendMsg
        od4.timeTrigger(FREQ,sendMsg);
    }
    return 0;
}

