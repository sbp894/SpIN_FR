function x_out= return_closest_max(x_in, y_in, x_test, x_search_spread, x_output_spread)

if ~exist('x_output_spread', 'var')
   x_output_spread=0; 
end

x_in= x_in(:);
y_in= y_in(:);

if ~exist('x_spread', 'var')
    x_search_spread=5;
end

x_out= nan(length(x_test)*(2*x_output_spread+1), 1);
for testVar=1:length(x_test)
    cur_test_ind= dsearchn(x_in, x_test(testVar));
    inds2use= max(1, cur_test_ind-x_search_spread):min(numel(x_in), cur_test_ind+x_search_spread);
    [~, indShift]= max(y_in(inds2use));
    indAbsolute= inds2use(1)-1+indShift;
    x_out(((testVar-1)*(2*x_output_spread+1)+1):(testVar*(2*x_output_spread+1)))= x_in((indAbsolute-x_output_spread) : (indAbsolute+x_output_spread));
end