#ifndef test_setup_h
#define test_setup_h

// Comment this line for profiling settings
//#define USE_MATLAB false

//
// Use these settings for comparing with MATLAB output
//
#ifdef USE_MATLAB

#define N_FRAMES 2
#define I_FRAME_FREQ 2
#define DUMP_TO_DEBUG true
#define image_name "blooper"

#endif

//
// Use these settings for profiling optimizaitons
//
#ifndef USE_MATLAB
#define N_FRAMES 10
#define I_FRAME_FREQ 1
#define DUMP_TO_DEBUG false
#define image_name "solar"
#endif

#endif