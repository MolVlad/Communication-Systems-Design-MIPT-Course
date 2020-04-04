clear all
close all

runParam = param();
Scenario2(runParam)
main(runParam)

if runParam.isMeasurement
    runParam.isPlotSpectrum = 0;
    
    clc
    
    runParam.angleRange = 90*[0 1];
    main(runParam)

    runParam.angleRange = 90*[-1 0];
    main(runParam)
end