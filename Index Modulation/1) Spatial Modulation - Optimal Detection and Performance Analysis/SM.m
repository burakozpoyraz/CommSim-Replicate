%% Function Information
% ARGUMENTS
% 1-) num_iterations: Number of iterations for Monte Carlo simulation
% 2-) Nt: Number of transmit antennas
% 3-) Nr: Number of receive antennas
% 4-) M: Modulation level
% 5-) SNRdB: Signal-to-noise ratio in dB scale
% 6-) mod_type: Constellation scheme ("PSK" or "QAM")
% 7-) ch_type: Channel type ("Conventional" or "Constrained")
% 8-) det_type: Detector type ("Optimal" or "Mesleh")

% OUTPUTS
% 1-) BER: Bit error rate
% 2-) num_bit_errors: Number of bit errors

%% Main Function
function [BER, num_bit_errors] = SM(num_iterations, Nt, Nr, M, SNRdB, mod_type, ch_type, det_type)
    % Parameters///////////////////////////////////////////////////////////
    % APM Symbol Set (Constellation)---------------------------------------
    if mod_type == "QAM"
        ss = qammod(0 : M-1, M, "Gray", "UnitAveragePower", true);
    elseif mod_type == "PSK"
        ss = pskmod(0 : M-1, M, 0, "Gray");
    else
        error("The 'mod_type' argument must be either 'PSK' or 'QAM'.")
    end
    % ---------------------------------------------------------------------

    ns = log2(Nt); % Number of spatial bits for antenna selection
    nd = log2(M); % Number of data bits for APM symbol selection
    n_tot = ns + nd; % Number of total bits transmitted in a time slot
    num_bits = num_iterations * n_tot; % Number of total bits transmitted 
                                       % in whole simulation

    SNR = 10^(SNRdB / 10); % Signal-to-noise ratio in linear scale
    % /////////////////////////////////////////////////////////////////////

    % Simulation///////////////////////////////////////////////////////////
    num_bit_errors = 0;
    for iter_index = 1 : num_iterations
        % Transmitter------------------------------------------------------
        data_bits = randi([0, 1], [1, nd]);
        s = ss(Bit2Dec(data_bits) + 1); % Transmitted APM symbol

        spatial_bits = randi([0, 1], [1, ns]);
        i = Bit2Dec(spatial_bits) + 1; % Activated antenna index

        x = zeros(Nt, 1);
        x(i) = s; % Transmitted signal array
        % -----------------------------------------------------------------

        % Channel & Noise--------------------------------------------------
        % Rayleigh fading channel with i.i.d. entries according to CN(0, 1)
        if ch_type == "Conventional"
            H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
        elseif ch_type == "Constrained"
            H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
            norm_array = vecnorm(H, 2, 1);
            H = H ./ norm_array;
        else
            error("The 'ch_type' argument must be either 'Conv' or 'Cons'.")
        end

        % AWGN with i.i.d. entries according to CN(0, 1)
        n = (randn(Nr, 1) + 1i * randn(Nr, 1)) / sqrt(2);
        % -----------------------------------------------------------------

        % Receiver---------------------------------------------------------
        y = sqrt(SNR) * H * x + n;
        [detected_data_bits, detected_spatial_bits] = ...
            Detector(ss, H, y, Nt, M, SNR, det_type);

        num_data_bit_errors = sum(xor(data_bits, detected_data_bits));
        num_spatial_bit_errors = sum(xor(spatial_bits, detected_spatial_bits));
        num_bit_errors = num_bit_errors + num_data_bit_errors + num_spatial_bit_errors;
        % -----------------------------------------------------------------
    end
    BER = num_bit_errors / num_bits;
    % /////////////////////////////////////////////////////////////////////
end

%% Inner Functions
% =========================================================================
% 1. Conversion from bit to decimal
%
% ARGUMENT
% - bit_array: Bit array to be converted to decimal value
%
% OUTPUT
% - decimal_value: Corresponding decimal value of the bit array
% =========================================================================
function decimal_value = Bit2Dec(bit_array)
    size_bit_array = size(bit_array);
    num_bits = size_bit_array(2);
    decimal_value = bit_array * (2.^((num_bits - 1) : -1 : 0))';
end
% =========================================================================


% =========================================================================
% 2. Detecting transmitted data and spatial bits
%
% ARGUMENTS
% 1-) ss: APM symbol set
% 2-) H: Rayleigh fading channel
% 3-) y: Received signal array
% 4-) Nt: Number of transmit antennas
% 5-) M: Modulation level
% 6-) SNR: Signal-to-noise ratio in linear scale
% 7-) det_type: Detector type ("Optimal" or "Mesleh")
%
% OUTPUTS
% 1-) detected_data_bits: Detected data bits corresponding to the
%     transmitted APM symbol
% 2-) detected_spatial_bits: Detected spatial bits corresponding to the
%     activated antenna index
% =========================================================================
function [detected_data_bits, detected_spatial_bits] = Detector(ss, H, y, Nt, M, SNR, det_type)
    nd = log2(M);
    ns = log2(Nt);

    if det_type == "Optimal"
        error_matrix = zeros(Nt * M, 3);
        comb_index = 1;
        for i = 1 : Nt
            hi = H(:, i);
            for k = 1 : M
                sk = ss(k);

                error_matrix(comb_index, 1) = i;
                error_matrix(comb_index, 2) = k;

                error_value = norm(y - sqrt(SNR) * hi * sk, "fro")^2;
                error_matrix(comb_index, 3) = error_value;
                comb_index = comb_index + 1;
            end
        end
        error_array = error_matrix(:, end);
        [~, min_error_index] = min(error_array);

        detected_i = error_matrix(min_error_index, 1);
        detected_spatial_bits = Dec2Bit(detected_i - 1, ns);

        detected_s_index = error_matrix(min_error_index, 2);
        detected_data_bits = Dec2Bit(detected_s_index - 1, nd);
    elseif det_type == "Mesleh"
        z_array = zeros(1, Nt);
        for i = 1 : Nt
            hi = H(:, i);
            zi = hi' * y / sqrt(SNR);
            z_array(i) = zi;
        end
        [~, detected_i] = max(z_array);
        detected_spatial_bits = Dec2Bit(detected_i - 1, ns);

        zi = z_array(detected_i);
        [~, detected_s_index] = min(abs(zi - ss).^2);
        detected_data_bits = Dec2Bit(detected_s_index - 1, nd);
    else
        error("The 'det_type' argument must be either 'Optimal' or 'Mesleh'.")
    end
end
% =========================================================================


%==========================================================================
% 3. Conversion from Decimal to Bit

% ARGUMENTS
% 1-) decimal: Decimal value to be converted to bit array
% 2-) n: Number of bits that the resulting bit array has

% OUTPUT
% - bit_array: Corresponding bit array of the decimal value
%==========================================================================
function bit_array = Dec2Bit(decimal, n)
    decreasing_pow_array = (n - 1 : -1 : 0);
    bit_array = zeros(length(decimal), n);
    for bit_index = 1 : n
        pow_val = decreasing_pow_array(bit_index);
        comparison = (decimal / 2^pow_val) >= 1;
        false_array = find(comparison == 0);
        true_array = find(comparison == 1);
        bit_array(false_array, bit_index) = 0;
        bit_array(true_array, bit_index) = 1;
        decimal(true_array) = decimal(true_array) - 2^pow_val;
    end
end
%==========================================================================
