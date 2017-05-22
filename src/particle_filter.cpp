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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;	//defining number of particles
	weights.resize(num_particles);
	particles.resize(num_particles);

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
		//cout << "particle: " << particles[i].x << endl;
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

		if (fabs(yaw_rate) < 0.00001) {
			particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
		}
		else {
			particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta = particles[i].theta + yaw_rate*delta_t;
		}

		//cout << "x " << particles[i].x <<endl;

		// creates a normal (Gaussian) distribution for x,y,theta
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

		//cout << "x noise " << particles[i].x << endl;

	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	int nearest_neighbor;
	double po_dist_past = -1.0; //initialize the previous calculated distance
	for (int i=0; i<observations.size(); i++) {    //loop thru predicted landmarks
		for (int j=0; j<predicted.size(); j++) {    //loop through landmark observations
			double po_dist = dist(observations[i].x , observations[i].y , predicted[j].x , predicted[j].y);
			if (po_dist_past < 0) {
				nearest_neighbor = predicted[j].id; //select po_dist as the nearest neighbor
				po_dist_past = po_dist;
			}
			else if (po_dist < po_dist_past) {
				nearest_neighbor = predicted[j].id; //select po_dist as the new nearest neighbor.
				po_dist_past = po_dist;
			}
				
		}
		//cout << "obs id before: " << observations[i].id << endl;
		observations[i].id = nearest_neighbor;	//assign the id of the nearest predicted landmark to the observation
		//cout << "obs id after: " << observations[i].id << endl;
		po_dist_past = -1.0;
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

	weights.clear();
	

	for (int i=0; i<num_particles; i++) {    //loop thru particles
		particles[i].weight = 1; //reinitialize at 1
		std::vector<LandmarkObs> tobservations;
		for (int j=0; j<observations.size(); j++) {    //loop thru observations
			double tobsx = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;	//transform into map coords
			double tobsy = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;
			tobservations.push_back(LandmarkObs{observations[j].id,tobsx,tobsy});
			//cout << "p x,y,th: " << particles[i].x <<" "<<particles[i].y<<" "<<particles[i].theta<<endl;
			//cout << "obs x: " << observations[j].x << endl;
		}
		//call data association to place a landmark ID with each observation.
		std::vector<LandmarkObs> pred_map_data;
		for (int k=0; k<map_landmarks.landmark_list.size(); k++) {
			double landmark_dist = dist(map_landmarks.landmark_list[k].x_f , map_landmarks.landmark_list[k].y_f , particles[i].x  , particles[i].y);
			if (landmark_dist <= sensor_range) {
				LandmarkObs mark_temp;
				mark_temp.id = map_landmarks.landmark_list[k].id_i; //check indices
				mark_temp.x = map_landmarks.landmark_list[k].x_f;
				mark_temp.y = map_landmarks.landmark_list[k].y_f;
				pred_map_data.push_back(mark_temp);
			}
		}
		
		dataAssociation(pred_map_data, tobservations); //associate a landmark with each observation
		//cout << "tobs: " << tobservations[1].id << " " << tobservations[1].x << " "<<tobservations[1].y << endl;
		//cout << "pred: " << map_landmarks.landmark_list[(tobservations[1].id -1)].id_i << " "<< map_landmarks.landmark_list[(tobservations[1].id -1)].x_f <<" "<< map_landmarks.landmark_list[(tobservations[1].id -1)].y_f << endl;
		//calculate probabilities using multivariate gaussian probability for each observation, multiply with the weight.
		float sigx = std_landmark[0];
		float sigy = std_landmark[1];
		for (int j=0; j<tobservations.size(); j++) {
			int mark_num = tobservations[j].id -1;
			//cout << "mark_num: " << mark_num << endl;
			//cout << "landmark: " << map_landmarks.landmark_list[mark_num].id_i <<endl;
			float x_sq = pow((tobservations[j].x - map_landmarks.landmark_list[mark_num].x_f),2);
			float y_sq = pow((tobservations[j].y - map_landmarks.landmark_list[mark_num].y_f),2);
			//cout << "x_sq= " << x_sq << " y_sq= " << y_sq << endl;
			float prb = (1/(2*3.1415*sigx*sigy))*exp(-(x_sq/(2*sigx*sigx)+y_sq/(2*sigy*sigy)));
			//cout << "prb= " << prb << endl;
			particles[i].weight *= prb;
		}
		//cout << particles[i].weight << endl; debugging purposes only
		weights.push_back(particles[i].weight);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles and replace with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::default_random_engine gen;
    std::vector<Particle> resampled_particles;
    resampled_particles.resize(num_particles);
    std::discrete_distribution<> d(weights.begin(), weights.end());

	for(int i=0; i<num_particles; ++i) {
        resampled_particles[i] = particles[d(gen)];
    }
    particles=resampled_particles;
	
	/* for debugging purposes only
	for (int i=0; i<num_particles; i++) {
		cout << "weight " << i+1 << ": " << particles[i].weight << endl;
	}
	*/

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
