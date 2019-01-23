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
 * This header file contains necessary data structures
 * for code transplant from Matlab to C++
 *
 * Yue Kang and Yushu Yu, 
 * Revere Lab, Chalmers / GU, 2018
*/

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cstdlib>
#include <ctime>
#include <cmath>

class FB_state
{
public:
    double xp_dot{0.0};  // longitudinal speed
    double yp_dot{0.0};  // lateral speed
    double psi_dot{0.0}; //
    double epsi{0.0};    //
    double ey{0.0};      // lateral position
    double s{0.0};       // longitudinal position
    double steer{0.0};
    double acc{0.0};
    FB_state() {}
    FB_state(double a, double b, double c, double d, double e, double f, double g, double h):
        xp_dot(a), yp_dot(b), psi_dot(c), epsi(d), ey(e), s(f), steer(g), acc(h){}
    FB_state operator +(const FB_state& other)
    {
        FB_state new_state(
            this->xp_dot + other.xp_dot,
            this->yp_dot + other.yp_dot,
            this->psi_dot + other.psi_dot,
            this->epsi + other.epsi,
            this->ey + other.ey,
            this->s + other.s,
            this->steer + other.steer,
            this->acc + other.acc
        );
        return new_state;
    }
    FB_state operator +(const Eigen::VectorXd& other)
    {
        FB_state new_state( 
            this->xp_dot + other(0),
            this->yp_dot + other(1),
            this->psi_dot + other(2),
            this->epsi + other(3),
            this->ey + other(4),
            this->s + other(5),
            this->steer + other(6),
            this->acc + other(7)
        );
        return new_state;
    }
    FB_state& operator +=(const FB_state& other)
    {
        this->xp_dot += other.xp_dot;
        this->yp_dot += other.yp_dot;
        this->psi_dot += other.psi_dot;
        this->epsi += other.epsi;
        this->ey += other.ey;
        this->s += other.s;
        this->steer += other.steer;
        this->acc += other.acc;
        return *this;
    }
    FB_state& operator +=(const Eigen::VectorXd& other)
    {
        this->xp_dot += other(0);
        this->yp_dot += other(1);
        this->psi_dot += other(2);
        this->epsi += other(3);
        this->ey += other(4);
        this->s += other(5);
        this->steer += other(6);
        this->acc += other(7);
        return *this;
    }
    void print()
    {
        std::cout << "[" << this->xp_dot << ", " << this->yp_dot << ", " << this->psi_dot << ", " << this->epsi;
        std::cout << ", " << this->ey << ", " << this->s << ", " << this->steer << ", " << this->acc << "]" << std::endl;
    }
};

class Obstacle
{
public:
    double pos_x{0.0};
    double pos_y{0.0};
    double vel_x{0.0};
    double vel_y{0.0};
    double acc_x{0.0};
    double acc_y{0.0};
    double radius{0.0};

    bool isConf(Obstacle that)
    {
        Eigen::Vector2d temp;
        temp << this->pos_x - that.pos_x, this->pos_y - that.pos_y;
        return (temp.norm() <= (this->radius + that.radius)) ? true : false;
    }

    bool isConf(std::vector<Obstacle> list)
    {
        if (0 == list.size()) return false;
        for (uint16_t i = 0; i < list.size(); i++)
        {
            if (this->isConf(list[i])) return true;
        }
    return false;
    }

    void print()
    {
        std::cout << "Obstacle coor: [" << pos_x << "; " << pos_y << "], vel: [" << vel_x << "; " << vel_y << "]" << std::endl;
        std::cout << "acc:  [" << acc_x << "; " << acc_y << "], radius: " << radius << std::endl;
    }
};

class Coefficient // the return value "out", line 463-480 in constraint_obs~.m
{
public:
    double norm_relpos{0.0};
    double h_angle_moving{0.0};
    Eigen::Vector2d A_n_angle_moving; 
    double b_n_angle_moving{0.0};    
    double h_angle_fix{0.0}; 
    Eigen::Vector2d A_n_angle_fix;
    double b_n_angle_fix{0.0};
    double h_dis{0.0};
    Eigen::Vector2d A_n_dis;
    double b_n_dis{0.0};
    bool alert{false};
    double h_sid_pos{0.0};
    Eigen::Vector2d A_n_side_pos;
    double b_n_side_pos{0.0};
    double h_sid_neg{0.0};
    Eigen::Vector2d A_n_side_neg;
    double b_n_side_neg{0.0};
    double radius{0.0};

    void print()
    {
        std::cout << "constraint(norm_relpos, angle):" << this->norm_relpos << ", " << this->h_angle_fix  << ", " << this->A_n_angle_fix(0) << ", "  << this->A_n_angle_fix(1)<< ", " << this->b_n_angle_fix << std::endl;
	std::cout << "constraint(distance): "  << this->h_dis<< ", "   <<  this->A_n_dis(0) << ", " << this->A_n_dis(1) << ", " << this->b_n_dis << std::endl;
	std::cout << "constraint(postive side): "  << this->h_sid_pos << ", "   <<  this->A_n_side_pos(0) << ", " << this->A_n_side_pos(1) << ", " << this->b_n_side_pos << std::endl;
	std::cout << "constraint(negative side): " << this->h_sid_neg  << ", "  <<  this->A_n_side_neg(0) << ", " << this->A_n_side_neg(1) << ", " << this->b_n_side_neg << std::endl;
    }
};

class Output_safety // the return value output_safety.out, line 355-359 in safety_cert~.m
{
public:
    Eigen::Vector2d x;
    Coefficient coef;
    double value_min{100000000.0};
    bool hasSolution{false}; // line 298 in safety_cert~.m, INVERSED BOOLEAN VALUE of "nosolution"
};

class Global_variables
{
public:
    int scale{0}; // mayby not useful
    Eigen::Vector2d u_global{};
    int scale_tracking{0};
    Eigen::Vector2d u_tracking_global{};
    int scale_record{0};

    std::vector<Eigen::Vector3d> trajd{}; // this contains all the 3 variables in the following line
    //Eigen::Vector3d tra_com_pre, tra_com_dot_pre, tra_com_ddot_pre;

    bool brake_flag{false}, brake_flag_pre{false};
    bool nosolution{false};

    // the following three are initialised by constraint_obstacles.m
    std::vector<Obstacle> traj_ob{};
    int no_ob{3};
    // vector<Eigen::Vector2d> pos_ob_array_pre;
    // vector<double> radius_pre;

    std::vector<bool> beta_2{};
    double dt{0.0};

    bool isVerbose{false};

    // // For data sample and visualisation
    // vector<double> t_ctrl;
    // vector<Eigen::VectorXd> u_ctrl;

    // Unlisted global variables
    bool dead{false};
    FB_state state_brakeini{};
    Global_variables()
    {
        scale = 0; 
        u_global << 0.0, 0.0;
        scale_tracking = 0;
        u_tracking_global << 0.0, 0.0;
        scale_record = 0;
        Eigen::Vector3d tra_com_pre, tra_com_dot_pre, tra_com_ddot_pre;
        tra_com_pre << 0.0, 0.0, 0.0;
        tra_com_dot_pre << 0.0, 0.0, 12.0; // Pay attention to the value here!
        tra_com_ddot_pre << 0.0, 0.0, 0.0;
        trajd.push_back(tra_com_pre);
        trajd.push_back(tra_com_dot_pre);
        trajd.push_back(tra_com_ddot_pre);
        brake_flag = false;
        brake_flag_pre = false;
        nosolution = false;
        dt = 0.001;
        for (int i = 0; i < 50; ++i)
        {
            beta_2.push_back(false);
        }
    }
    void generate_init_ob()
    {
        using namespace std;
        srand((int)time(NULL));
        for (uint16_t i = 0; i < no_ob; i++)
        {
            Obstacle curr;
            do
            {
                curr.radius = 1.5 + 0.5 * (double)rand() / RAND_MAX;
                curr.pos_x = 16 + 80 * (double)rand() / RAND_MAX;
                curr.pos_y = -1.2 + 2.4 * (double)rand() / RAND_MAX;

		//tune 20190110:
		/*curr.radius = 2.1; 
		curr.pos_x = 86.11;
		curr.pos_y = 1.342;*/

		//curr.radius = 2.1; 
		//curr.pos_x = 90.11;
		//curr.pos_y = 0.342;
            }
            while (curr.isConf(traj_ob));
            traj_ob.push_back(curr);
        }
    }

    void ob_traj(bool isDynamic) // update position of the obstacles
    {
        for (uint16_t i = 0; i < no_ob; i++)
        {
            double v = isDynamic ? std::sqrt(dt) * 2 : 0.0;

            traj_ob[i].pos_x += v * dt;
            traj_ob[i].vel_x = v;
        }
    }

    void traj_gen(FB_state y) // update trajd
    {
        double t = ((double)clock())/CLOCKS_PER_SEC;
        double ks;
        if (nosolution)
        {
            ks = 0.1;
            nosolution = false;
        }
        else
        {
            ks = 0.0;
        }
        Eigen::Vector3d err;
        err << y.epsi, y.ey, y.s;
        err -= trajd[0]; // tra_com_pre
        double virtual_time = std::exp(-ks * err.squaredNorm());
        double v = 20 * virtual_time;

	v = 12; //constant 
        //trajd[0](2) = v * t; // tra_com_pre
        trajd[1](2) = v; // tra_com_dot_pre
	trajd[0](1) = 0; // 
        trajd[1](1) = 0; // tra_com_dot_pre
    }
};


std::vector<Coefficient> constraint_obstacles_dynamics_complex(FB_state, Global_variables &);

Output_safety safety_certificate_complex(FB_state, Global_variables &);
