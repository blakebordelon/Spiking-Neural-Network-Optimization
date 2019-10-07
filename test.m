% plot the functions f(t) = (3/2)e-t-(1/2)e-3t and g(t) = e-2t+2te-2t for 0 ? t ? 5 on the same plot.
% 101 points will be plotted to get a smooth curve. deltat wil be the spacing between points.
	deltat = 5/100;
% store the time and function values in 3 vectors of length 200;
	for i  = 1:200
		t(i) = (i-1) * deltat;
		foft(i) = 3/2 *exp(-t(i))-1/2*exp(-3*t(i));
		goft(i) = exp(-2 * t(i)) + 2*t(i)*exp(-2*t(i));
	end
% plot and label the results
	plot(t, foft, 'r*', t, goft, '.b');
	grid;
	legend('f(t) = (3/2)e-t-(1/2)e-3t','g(t) = e-2t+2te-2t');
	xlim([0 10]);
	ylim([-1 2]);
	xlabel('time');
	ylabel('function');
