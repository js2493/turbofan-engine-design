function massflowrate = mdot(maxTemp,compRatio)
massflowrate=2.5*2*231.8*(maxTemp/216.65)/sqrt(compRatio);
end