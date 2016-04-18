struct timer {
	double tot_diff;
	time_t begin_time; // double
	time_t cur_time; // double
	double sec;
	int min;
	int hr;
	int time_depth; // {1,2,3} = only up to {sec, min, hr}
};

// Timer functions //
void start_timer(timer &my_timer, bool restart = 1);
void current_time(timer &my_timer);
void inc_time_diff(timer &my_timer);
double display_total_time_diff(timer &my_timer, int depth = 1);
double display_elapse_time(timer &my_timer, int depth = 1);
double display_remaining_time(timer &my_timer, double cur_iter,
                              double tot_iter,  int depth = 1);
