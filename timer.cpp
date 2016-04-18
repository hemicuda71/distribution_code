#include <iostream>
#include <iomanip>
//#include <mpi.h>
#include "timer.h"

using namespace std;

// Timer function definitions //
void start_timer(timer &my_timer, bool restart) {
	if(restart)
		my_timer.tot_diff = 0.0;
	//my_timer.begin_time = MPI_Wtime();
	time(&my_timer.begin_time);
}

void current_time(timer &my_timer) {
	//my_timer.cur_time = MPI_Wtime();
	time(&my_timer.cur_time);
}

void inc_time_diff(timer &my_timer) {
	//my_timer.tot_diff += my_timer.cur_time - my_timer.begin_time;
	my_timer.tot_diff += difftime(my_timer.cur_time, my_timer.begin_time);
}

double display_total_time_diff(timer &my_timer, int depth) {
	my_timer.time_depth = depth;
	//my_timer.sec = difftime(my_timer.cur_time, my_timer.begin_time);
	my_timer.sec = my_timer.tot_diff;
	my_timer.hr = my_timer.sec / 3600;
	my_timer.sec -= 3600 * my_timer.hr;
	my_timer.min = my_timer.sec / 60;
	my_timer.sec -= 60 * my_timer.min;
	cout << "Total time: ";
	switch (my_timer.time_depth) {
	case 3:
		cout << my_timer.hr << "hr  ";
	case 2:
		cout << my_timer.min << "min  ";
	case 1:
		cout << my_timer.sec << "sec\n";
	}
	return my_timer.tot_diff;
}

double display_elapse_time(timer &my_timer, int depth) {
	my_timer.time_depth = depth;
	my_timer.sec = difftime(my_timer.cur_time, my_timer.begin_time);
	//my_timer.sec = my_timer.cur_time - my_timer.begin_time;
	my_timer.hr = my_timer.sec / 3600;
	my_timer.sec -= 3600 * my_timer.hr;
	my_timer.min = my_timer.sec / 60;
	my_timer.sec -= 60 * my_timer.min;
	cout << "Time elapsed: ";
	switch (my_timer.time_depth) {
	case 3:
		cout << my_timer.hr << "hr  ";
	case 2:
		cout << my_timer.min << "min  ";
	case 1:
		cout << my_timer.sec << "sec\n";
	}
	return my_timer.cur_time - my_timer.begin_time;
}

double display_remaining_time(timer &my_timer, double cur_iter, double tot_iter, int depth)
{
	my_timer.time_depth = depth;
	my_timer.sec = difftime(my_timer.cur_time, my_timer.begin_time);
	//my_timer.sec = my_timer.cur_time - my_timer.begin_time;
	my_timer.sec *= (tot_iter / cur_iter - 1.0);
	my_timer.hr = my_timer.sec / 3600;
	my_timer.sec -= 3600 * my_timer.hr;
	my_timer.min = my_timer.sec / 60;
	my_timer.sec -= 60 * my_timer.min;
	cout << "Time remaining: ";
	switch (my_timer.time_depth) {
	case 3:
		cout << my_timer.hr << "hr  ";
	case 2:
		cout << my_timer.min << "min  ";
	case 1:
		cout << setprecision(4) << my_timer.sec << "sec\n";
	}
	return my_timer.cur_time - my_timer.begin_time;
}
