/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */


#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include "particle_filter.h"
using namespace std;

void ParticleFilter::cvrtUpdateOnParticle(vector<Particle> &Particles,double delta_t,double std_pos[], double velocity, double yaw_rate)
{
}

std::ostream& dump(std::ostream &o, const LandmarkObs& p)
{
    return o << "x: " << p.x << "\ty: " << p.y << "\tid" <<p.id << std::endl;
}


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 1;
    normal_distribution<double> N_x_init(0, std[0]);
    normal_distribution<double> N_y_init(0, std[1]);
    normal_distribution<double> N_theta_init(0, std[2]);
    for(auto i=0; i<num_particles;i++ ) {
        Particle p;
        double n_x = N_x_init(gen);
        double n_y = N_y_init(gen);
        double n_theta = N_theta_init(gen);
        p.x = x + n_x;
        p.y = y + n_y;
        p.theta = theta + n_theta;
        particles.push_back(p);
    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    normal_distribution<double> N_x_init(0, std_pos[0]);
    normal_distribution<double> N_y_init(0, std_pos[1]);
    normal_distribution<double> N_theta_init(0, std_pos[2]);
    for(auto &particle : particles) {
        double n_x = N_x_init(gen);
        double n_y = N_y_init(gen);
        double n_theta = N_theta_init(gen);
        if (fabs(yaw_rate) > 0.001) {

           cout <<"Test " <<  particle.theta * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta))<<"\n";

            particle.x = particle.x + velocity / yaw_rate * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta)) + n_x;
            particle.y = particle.y + velocity / yaw_rate * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t)) + n_y;
            particle.theta = particle.theta + yaw_rate * delta_t + n_theta;


           // cout << "This is executing";
        } else {
            particle.x = particle.x + velocity * delta_t * cos(particle.theta) + n_x;
            particle.y = particle.y + velocity * delta_t * sin(particle.theta) + n_y;
            particle.theta = particle.theta + n_theta;
            cout << "That is executing ";
        }
    }

      cout << "particle x : " << particles[0].x << "\t";
      cout << "particle y : " <<particles[0].y;
    // cout<< "\n";

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
       double min_distance, dist, dx, dy;

        for(auto const p :predicted)
        {
             min_distance = INFINITY;
            for(auto &o: observations)
            {

                dx = (p.x - o.x);
                dy = (p.y - o.y);
                dist = dx*dx + dy*dy;
                if( dist < min_distance ) {
                    o.id = p.id;
                    min_distance = dist;
                }
            }
        }
}

double multivariate_normal_distribution(LandmarkObs predicted, LandmarkObs observations,double std_landmark[]) {
    /* returns the probability of in mutli-variate gaussian distribution
     * predicted contains the mean and landmark is the value of x,y for which we need to find the probability
     */

    double coeff = (2.0*M_PI*std_landmark[0]*std_landmark[1]);
    double sigma_square_x = std_landmark[0] * std_landmark[0];
    double sigma_square_y = std_landmark[1] * std_landmark[1];
    double dx = (observations.x - predicted.x);
    double dy = (observations.y - predicted.y);

    double x_diff = dx * dx / (2 * sigma_square_x);
    double y_diff = dy * dy / (2 * sigma_square_y);


   double prob_x_y =  exp(-(x_diff+y_diff)) / coeff;
   //  cout << "Prob_xy : " << prob_x_y <<"\n";

    return prob_x_y;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html


    // if the observation is above the sensor range, ignore it
    weights.clear();
    for(auto &p : particles) {

        std::vector<LandmarkObs> landmark;
        for( auto const lm: map_landmarks.landmark_list)
        {
            auto dis_x = (lm.x_f - p.x) * (lm.x_f - p.x);
            auto dis_y = (lm.y_f - p.y) * (lm.y_f - p.y);
            if(sqrt(dis_x + dis_y) < sensor_range)
            {
                LandmarkObs l;
                l.x = lm.x_f;
                l.y  = lm.y_f;
                l.id = lm.id_i;
                cout<<"landmark";
                landmark.push_back(l);
                dump(cout, l);
            }

        }
        std::vector<LandmarkObs> obs_particle;

        // Transforming from particle coordinate to global coordinate
        for (auto const ob: observations) {
            LandmarkObs l_transformed;
            // to the observation of the particles

            l_transformed.x = ob.x * cos(p.theta) - ob.y * sin(p.theta) + p.x;
            l_transformed.y = ob.x * sin(p.theta) + ob.y * cos(p.theta) + p.y;
            l_transformed.id = ob.id;

            dump(cout,l_transformed);
            obs_particle.push_back(l_transformed);
        }
        // finds the nearest landmark associated with observation by the particle
        dataAssociation(landmark,obs_particle);
      //  std::cout<<"Car observation Size :"<< car_observation.size() <<"\n";
       // std::cout<<"Paticle observation size :"<< obs_particle.size() << "\n";
//        for (auto &ob: obs_particle) {
//            cout <<"particle observation id : "  << ob.id << " \n";
//        }
        double  weight = 1.0;
        for(int i=0; i < obs_particle.size(); i++) {

            for (auto const lm: landmark) {
                //  check for associated landmark, since multiple landmark can be associated with same landmark
                if (obs_particle[i].id == lm.id) {

                    weight = weight * multivariate_normal_distribution(lm, obs_particle[i], std_landmark);
                }
               // cout <<"obs_particle.id : " << obs_particle[i].id<<"\t";
            }
          //  cout << "\n";
        }
        cout<< "weight : " << weight << "\n";
        weights.push_back(weight);
        p.weight = weight;

    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::discrete_distribution<> d(weights.begin(),weights.end());
    std::vector<Particle> gen_particles;
    for(int i=0; i < num_particles; i++)
    {
        int ran = d(gen);
        //cout << "randome numbers: " << ran << "\t";
        gen_particles.push_back(move(particles[ran]));
    }
    cout<<"\n";
    particles = move(gen_particles);
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
