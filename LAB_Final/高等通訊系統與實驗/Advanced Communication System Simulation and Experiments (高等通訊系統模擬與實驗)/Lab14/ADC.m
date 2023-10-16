function ADC_output = ADC(value,NOB,DR)
    NQL = 2^NOB;
    step = 2*DR/NQL;
    quantize = [-DR:step:DR];
    quantize = quantize(1:end-1);
    distance = abs(quantize-value);
    [minimum,index] = min(distance);
    ADC_output = quantize(index);
end
