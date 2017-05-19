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

#define _USE_MATH_DEFINES

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1000;	//defining number of particles

	default_random_engine gen;
	// creates a normal (Gaussian) distribution for x,y,theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	//particles.resize(num_particles,0);
	for (int i=0; i<num_particles; i++) {
		particles[i].id = i+1;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		particles[i].weight = 1.0;
		weights.push_back(particles[i].weight);
	}

	is_initialized = true;	//after setting these, set initialized to true
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	for (int i=0; i<num_particles; i++) {
		particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
		particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
		particles[i].theta = particles[i].theta + yaw_rate*delta_t;

		// creates a normal (Gaussian) distribution for x,y,theta
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

		/* maybe need this:
		normal_distribution<double> noise_x(0, std_pos[0]);
		normal_distribution<double> noise_y(0, std_pos[1]);
		normal_distribution<double> noise_theta(0, std_pos[2]);

		particles[i].x = particles[i].x + noise_x(gen);
		particles[i].y = particles[i].y + noise_y(gen);
		particles[i].theta = particles[i].theta + noise_theta(gen);
		*/
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	float nearest_neighbor;
	float po_dist_past = -1.0; //initialize the previous calculated distance
	for (int i=0; i<observations.size(); i++) {    //loop thru predicted landmarks
		for (int j=0; j<predicted.size(); j++) {    //loop through landmark observations
			float po_dist = sqrt(pow((observations[i].x - predicted[j].x),2) + pow((observations[i].y - predicted[j].y),2));
			if (po_dist_past < 0) {
				nearest_neighbor = predicted[j].id; //select po_dist as the nearest neighbor
			}
			else if (po_dist < po_dist_past) {
				nearest_neighbor = predicted[j].id; //select po_dist as the new nearest neighbor.
			}
				
		}
		observations[i].id = nearest_neighbor;	//assign the id of the nearest predicted landmark to the observation
	}
	//so now the loop is done. each osbservation will have a nearest neighbor landmark, assigned in its id
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
	
	for (int i=0; i<num_particles; i++) {    //loop thru particles
		for (int j=0; j<observations.size(); j++) {    //loop thru observations
			observations[j].x = cos(particles[i].theta)*observations[j].x + sin(particles[i].theta)*observations[j].y + particles[i].x;	//transform into map coords
			observations[j].y = sin(particles[i].theta)*observations[j].x - cos(particles[i].theta)*observations[j].y + particles[i].y;
		}
		//call data association to place a landmark ID with each observation.
		std::vector<LandmarkObs> pred_map_data;
		for (int k=0; k<map_landmarks.landmark_list.size(); k++) {
			float landmark_dist = sqrt(pow((map_landmarks.landmark_list[k].x_f - particles[i].x),2) + pow((map_landmarks.landmark_list[k].y_f - particles[i].y),2));
			if (landmark_dist < sensor_range) {
				LandmarkObs mark_temp;
				mark_temp.id = map_landmarks.landmark_list[k].id_i; //check indices
				mark_temp.x = map_landmarks.landmark_list[k].x_f;
				mark_temp.y = map_landmarks.landmark_list[k].y_f;
				pred_map_data.push_back(mark_temp);
			}
		}
		dataAssociation(pred_map_data, observations); //associate a landmark with each observation
		//calculate probabilities using multivariate gaussian probability for each observation, multiply with the weight.
		float sigx = std_landmark[0];
		float sigy = std_landmark[1];
		for (int j=0; j<observations.size(); j++) {
			int mark_num = observations[j].id -1;
			float prb = (1/(2*M_PI*sigx*sigy))*exp(-(pow((observations[j].x - map_landmarks.landmark_list[mark_num].x_f),2)/(2*sigx*sigx)
													+pow((observations[j].y - map_landmarks.landmark_list[mark_num].y_f),2)/(2*sigy*sigy)));
			particles[i].weight *= prb;
		}

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles and replace with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;	//set up random engine
	discrete_distribution<> d (weights); //create a weighted distribution of indices for which to resample particles from
	
	for (int i=0; i<num_particles; i++) {
		particles[i] = particles[d(gen)];
		weights[i] = particles[i].weight;
	}

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
