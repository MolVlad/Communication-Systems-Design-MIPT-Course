a = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];

data = [a a a a a a];

data = awgn(data, 10);

[r, lag] = xcorr(a, data);

stem(lag, r);