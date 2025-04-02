#include <avr/pgmspace.h>

// LSTM weights (int8, stored in PROGMEM)
const int8_t lstm_weight_ih[] PROGMEM = { /*...from lstm_weight_ih_int8.csv...*/ };
const int8_t lstm_weight_hh[] PROGMEM = { /*...from lstm_weight_hh_int8.csv...*/ };
const int8_t lstm_bias[] PROGMEM = { /*...from lstm_bias_int8.csv...*/ };

// Dense layers (int8)
const int8_t dense1_weight[] PROGMEM = { /*...*/ };
const int8_t dense1_bias[] PROGMEM = { /*...*/ };
// ...