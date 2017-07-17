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

	default_random_engine gen;

	num_particles = 100;

	normal_distribution<double> r_x(x, std[0]);
	normal_distribution<double> r_y(y, std[1]);
	normal_distribution<double> r_theta(theta, std[2]);


	for( int i = 0; i < num_particles; i++){
		Particle part;
		part.id = i;
		part.x = r_x(gen);
		part.y = r_y(gen);
		part.theta = r_theta(gen);
		part.weight = 1.0;

		weights.push_back(part.weight);
		particles.push_back(part);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;


	for( int i = 0; i < num_particles; i++){
		double pred_x;
		double pred_y;
		double pred_theta;

		if( yaw_rate == 0){
			pred_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			pred_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			pred_theta = particles[i].theta;
		}
		else{
			pred_x = particles[i].x + velocity/yaw_rate * ( sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta) );
			pred_y = particles[i].y + velocity/yaw_rate * ( cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t) );
			pred_theta = particles[i].theta + yaw_rate*delta_t;
		}

		// add jitter
		normal_distribution<double> r_x(pred_x, std_pos[0]);
		normal_distribution<double> r_y(pred_y, std_pos[1]);
		normal_distribution<double> r_theta(pred_theta, std_pos[2]);

		particles[i].x = r_x(gen);
		particles[i].y = r_y(gen);
		particles[i].theta = r_theta(gen);

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for( int i = 0; i < observations.size(); i++){
		// get first dist
		double x_dist = observations[i].x - predicted[0].x;
		double y_dist = observations[i].y - predicted[0].y;
		double curr_min_dist = sqrt ((x_dist*x_dist)+(y_dist*y_dist));
		double dist = curr_min_dist;
		// cout << curr_min_dist << endl;
		observations[i].id = predicted[0].id;
		for( int j = 1; j < predicted.size(); j++){
			x_dist = observations[i].x - predicted[j].x;
			y_dist = observations[i].y - predicted[j].y;
			dist = sqrt ((x_dist*x_dist)+(y_dist*y_dist));
			// cout << dist << endl;
			if( dist < curr_min_dist){
				// cout << "SMALLER DIST" << endl;
				// cout << dist << endl;
				curr_min_dist = dist;
				observations[i].id = j;
			}
		}
	}

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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	for (int i = 0; i < num_particles; i++){

		// transform observations

		vector<LandmarkObs> transformed_obs;
		LandmarkObs obs;
		for( int j = 0; j < observations.size(); j++){
			LandmarkObs temp_trans_obs;
			temp_trans_obs.x = particles[i].x + (observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta));
			temp_trans_obs.y = particles[i].y + (observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta));

			transformed_obs.push_back(temp_trans_obs);
		}

		vector<LandmarkObs> predicted;

		for( int j = 0; j < map_landmarks.landmark_list.size(); j++){
			predicted.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f});
		}

		dataAssociation(predicted, transformed_obs);

		particles[i].weight = 1.0;
		cout << transformed_obs.size() << endl;
		for( int j = 0; j < transformed_obs.size(); j++){
			double obs_x = transformed_obs[j].x;
			double obs_y = transformed_obs[j].y;
			double mu_x = map_landmarks.landmark_list[transformed_obs[j].id].x_f;
			double mu_y = map_landmarks.landmark_list[transformed_obs[j].id].y_f;
			double s_x = std_landmark[0];
			double s_y = std_landmark[1];
			long double weight_multiplier = 1/(2*M_PI*s_x*s_y)*exp(-((((obs_x-mu_x)*(obs_x-mu_x))/(2*s_x*s_x))+(((obs_y-mu_y)*(obs_y-mu_y))/(2*s_y*s_y))));
		
			// cout << weight_multiplier << endl;
			if( weight_multiplier > 0 ){
				particles[i].weight*=weight_multiplier;
			}
			

		}
		weights[i] = particles[i].weight;

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resampled_parts;

	for( int i = 0; i < num_particles; i++){
		resampled_parts.push_back(particles[distribution(gen)]);
	}

	particles = resampled_parts;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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
