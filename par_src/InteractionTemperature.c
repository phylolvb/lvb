/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Chang Sik Kim,
Maximilian Strobl and Martyn Winn
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/* ========== InteractionTemperature.c - MPI process management ========== */

#include "InteractionTemperature.h"

	/* release memory to the structrure */
	IterationTemperature *get_alloc_main_calc_iterations(void){
		IterationTemperature *p_temp;
		p_temp = (IterationTemperature *) alloc(sizeof(IterationTemperature), "alloc mem for IterationTemperature");
		p_temp->p_temperature = NULL;
		p_temp->p_next_iteration = NULL;
		p_temp->n_iteration = -1;
		return p_temp;
	}

	IndividualTemperature *get_alloc_temperature(SendInfoToMaster *p_info_temp, int n_try_process){
		IndividualTemperature *p_temp;
		p_temp = (IndividualTemperature *) alloc(sizeof(IndividualTemperature), "alloc mem for IndividualTemperature");
		p_temp->p_next_temperature = NULL;
		p_temp->d_temperature = p_info_temp->temperature;
		p_temp->n_seed = p_info_temp->n_seed;
		p_temp->l_length = p_info_temp->l_length;
		p_temp->n_try_process = n_try_process;
		return p_temp;
	}


	void add_temperature(IndividualTemperature *p_temperature, SendInfoToMaster *p_info_temp, int n_try_process){
		if (p_temperature->p_next_temperature == NULL){
			p_temperature->p_next_temperature = get_alloc_temperature(p_info_temp, n_try_process);
			return;
		}
		if (p_temperature->n_try_process == n_try_process) return; /* this doesn't occur */
		add_temperature(p_temperature->p_next_temperature, p_info_temp, n_try_process);
	}

	/* add temperature and info to structure */
	void add_temperature_cal_iterations(IterationTemperature *p_data, SendInfoToMaster *p_info_temp, int n_try_process){
		/* first added data */
		if (p_data->n_iteration == -1 || p_data->n_iteration == p_info_temp->n_iterations){
			p_data->n_iteration = p_info_temp->n_iterations;

			if (p_data->p_temperature == NULL) p_data->p_temperature = get_alloc_temperature(p_info_temp, n_try_process);
			else add_temperature(p_data->p_temperature, p_info_temp, n_try_process);
		}
		else{	/* next iteration temperature*/
			if (p_data->p_next_iteration == NULL) p_data->p_next_iteration = get_alloc_main_calc_iterations();
			add_temperature_cal_iterations(p_data->p_next_iteration, p_info_temp, n_try_process);
		}
	}



	double get_threshold_temperature(IndividualTemperature *p_temperature, int n_max_number_process){

		int n_count = 0;
		double d_temp = 0.0, d_avg, d_std;

		/* calc average */
		for (IndividualTemperature *p_temp = p_temperature; p_temp != NULL; p_temp = p_temp->p_next_temperature){
			n_count += 1;
			d_temp += p_temp->d_temperature;
		}
		if (n_count < n_max_number_process || n_count < 1) return 0.0;
		d_avg = d_temp / (double) n_count;

		/* calc std */
		d_temp = 0.0;
		for (IndividualTemperature *p_temp = p_temperature; p_temp != NULL; p_temp = p_temp->p_next_temperature){
			d_temp += pow_wrapper(d_avg - p_temp->d_temperature, 2);
		}
		d_std = sqrt(d_temp/(double) n_count);
		return d_avg + (d_std * (double) CALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS);
	}

	double get_threshold_length(IndividualTemperature *p_temperature, int n_max_number_process){

		int n_count = 0;
		long l_len = 0;
		double l_avg, l_std, l_ss;

		/* calc average */
		for (IndividualTemperature *p_temp = p_temperature; p_temp != NULL; p_temp = p_temp->p_next_temperature){
			n_count += 1;
			l_len += p_temp->l_length;
		}
		if (n_count < n_max_number_process || n_count < 1) return 0.0;
		l_avg = (double) l_len / (double) n_count;

		/* calc std */
		l_ss = 0.0;
		for (IndividualTemperature *p_temp = p_temperature; p_temp != NULL; p_temp = p_temp->p_next_temperature){
			//printf("\nerr:%f, result: %f\n", l_avg - (double)(p_temp->l_length), pow_wrapper(fabs(l_avg - (double)(p_temp->l_length)), 2));
			l_ss += pow_wrapper(fabs(l_avg - (double)(p_temp->l_length)), 2);
		}
		l_std = sqrt(l_ss / (double) n_count);
		return l_avg + (l_std * (double) CALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS);
	}

	
	Lvb_bool is_possible_to_continue(IterationTemperature *p_data, double d_temperature, int n_iteration,
				int l_tree_length, int n_max_number_process, int n_count_call){

		if (p_data == NULL) return LVB_TRUE;	/* default is always true */
		if (p_data->n_iteration == n_iteration) {	/* is this one to get the std + average */
			if (n_count_call < CALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS) return LVB_TRUE;

			double d_threshold_value = get_threshold_length(p_data->p_temperature, n_max_number_process);
	//		printf("Iteration:%d     threshold:%.8f    temperature:%.8f\n", n_iteration, d_threshold_value, d_temperature);
			if (d_threshold_value == 0)  return LVB_TRUE;
			if ((double) l_tree_length > d_threshold_value) return LVB_FALSE;
			return LVB_TRUE;
		}
		return is_possible_to_continue(p_data->p_next_iteration, d_temperature, n_iteration, l_tree_length, n_max_number_process, n_count_call + 1);
	}


	/* release memory to temperature structrure */
	void free_next_temperature(IndividualTemperature *p_temperature){
		if (p_temperature == NULL) return;
		free_next_temperature(p_temperature->p_next_temperature);
		p_temperature->p_next_temperature = NULL;
		free(p_temperature);
	}

	/* release memory to the structrure */
	void release_main_calc_iterations(IterationTemperature *p_data){

		if (p_data == NULL) return;
		release_main_calc_iterations(p_data->p_next_iteration);
		p_data->p_next_iteration = NULL;
		free_next_temperature(p_data->p_temperature);
		free(p_data);
	}
