/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 1000;
	std::default_random_engine generator;
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);
    for(int i=0;i<num_particles;i++){
        Particle particle;
        particle.id = i;
        particle.x = dist_x(generator);
        particle.y = dist_y(generator);
        particle.theta = dist_theta(generator);
        particle.weight = 1;
        particles.push_back(particle);
    }
    is_initialized = true;
    weights.reserve(num_particles);
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
 	std::default_random_engine generator;
    std::normal_distribution<double> dist_x(0, std_pos[0]);
    std::normal_distribution<double> dist_y(0, std_pos[1]);
    std::normal_distribution<double> dist_theta(0, std_pos[2]);
    for(int i=0;i<num_particles;i++){
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;
        if(fabs(theta)<0.001){
            particles[i].x = x + velocity*cos(theta)*delta_t ;
            particles[i].y = y + velocity*sin(theta)*delta_t ;
            //particles[i].theta = theta + yaw_rate*delta_t ;
        }else{
            particles[i].x = x + velocity/yaw_rate*(sin(theta+yaw_rate*delta_t)-sin(theta));
            particles[i].y = y + velocity/yaw_rate*(cos(theta) - cos(theta+yaw_rate*delta_t));
            particles[i].theta = theta + yaw_rate*delta_t;
        }
        particles[i].x += dist_x(generator);
        particles[i].y += dist_y(generator);
        particles[i].theta += dist_theta(generator);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
    for(unsigned int i=0;i<observations.size();i++){
        double x_obs = observations[i].x;
        double y_obs = observations[i].y;
        double min = numeric_limits<double>::max();
        for(unsigned int j=0;j<predicted.size();j++){
            double x_pred = predicted[j].x;
            double y_pred = predicted[j].y;
            int id = predicted[j].id;
            double d = dist(x_obs,y_obs,x_pred,y_pred);
            if(d<min){
                min = d;
                observations[i].id = id;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    double sum = 0.0;
    for(int k=0;k<num_particles;k++){
        double x_particle_map = particles[k].x;
        double y_particle_map = particles[k].y;
        double theta_particle = particles[k].theta;
        std::vector<LandmarkObs> predicted;
        for(unsigned int j=0;j<map_landmarks.landmark_list.size();j++){
            double x_landmark_map = map_landmarks.landmark_list[j].x_f;
            double y_landmark_map = map_landmarks.landmark_list[j].y_f;
            int id = map_landmarks.landmark_list[j].id_i;
            double dx = x_landmark_map - x_particle_map;
            double dy = y_landmark_map - y_particle_map;
            double landmark_range = sqrt(dx*dx+dy*dy);
            if( landmark_range <= sensor_range){
                LandmarkObs mark;
                mark.id = id;
                mark.x = x_landmark_map;
                mark.y = y_landmark_map;
                predicted.push_back(mark);
            }
        }
        std::vector<LandmarkObs> obs;
        for(int i=0;i<observations.size();i++){
            double x_obs_car = observations[i].x;
            double y_obs_car = observations[i].y;
            double x_obs_map = x_particle_map + x_obs_car*cos(theta_particle) - y_obs_car*sin(theta_particle);
            double y_obs_map = y_particle_map + x_obs_car*sin(theta_particle) + y_obs_car*cos(theta_particle);
            LandmarkObs ob;
            ob.x = x_obs_map;
            ob.y = y_obs_map;
            ob.id = observations[i].id;
            obs.push_back(ob);
        }

        dataAssociation(predicted,obs);
        double weight = 1.0;
        std::vector<int> associations;
        std::vector<double> sense_x;
        std::vector<double> sense_y;
        for(unsigned int i=0;i<obs.size();i++){
            double x_obs_map = obs[i].x;
            double y_obs_map = obs[i].y;

            int id = obs[i].id;
            for(int j=0;j<predicted.size();j++){
               if(id == predicted[j].id){
                   double x_landmark_map = predicted[j].x;
                   double y_landmark_map = predicted[j].y;
                   double deno = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
                   double nume = exp(
                         -(pow((x_obs_map-x_landmark_map),2)/(2*std_landmark[0]*std_landmark[0])+
                          pow((y_obs_map-y_landmark_map),2)/(2*std_landmark[1]*std_landmark[1]))
                          );
                   double weight_i = nume*deno;
                   weight *= weight_i;
                   associations.push_back(id);
                   sense_x.push_back(x_obs_map);
                   sense_y.push_back(y_obs_map);
                   //sense_x.push_back(x_landmark_map);
                   //sense_y.push_back(y_landmark_map);
                   break;
               }
            }
        }
        particles[k].weight = weight;
        //sum += particles[k].weight;
        SetAssociations(particles[k],associations,sense_x,sense_y);
    }
    for(int i=0;i<num_particles;i++){
       // particles[i].weight /= sum;
        weights[i]=particles[i].weight;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
 	std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,num_particles-1);
    int index = distribution(generator);
    double beta = 0.0;
    std::vector<double>::iterator max = std::max_element(
            weights.begin(),weights.end());
    double mw = *max;

    std::uniform_real_distribution<double> dist(0,1.0);
    std::vector<Particle> vparticle(particles);
    for(int i=0;i<num_particles;i++){
        beta += dist(generator)*2*mw;
        while(beta > weights[index]){
            beta -= weights[index];
            index = (index+1)%num_particles;
        }
        particles[i] = vparticle[index];
    }
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
