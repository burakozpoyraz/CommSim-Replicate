clear all;
clc;

Nt = 4;
Nr = 4;
M = 2;
SNRdB = 12;
mod_type = "PSK";
ch_type = "Cons";
det_type = "Mesleh";

ns = log2(Nt);
nd = log2(M);
n_tot = ns + nd;
num_min_bit_errors = 1100;
expected_BER = 0.015;
num_iterations = num_min_bit_errors / (n_tot * expected_BER);

[BER, num_bit_errors] = ...
    SM(num_iterations, Nt, Nr, M, SNRdB, mod_type, ch_type, det_type);