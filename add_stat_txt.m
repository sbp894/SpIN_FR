function add_stat_txt(mdl, params)
x_txt_val= params.x_txt_val;
y_txt_val= params.y_txt_val;
y_txt_gap= params.y_txt_gap;
fSize= params.fSize;
pValThresh= params.pValThresh;


if ~isempty(params.title)
    text(x_txt_val,y_txt_val,sprintf('%s', params.title), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex', 'FontWeight', 'bold');
    text(x_txt_val,y_txt_val-y_txt_gap,sprintf('%s', '------'), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex', 'FontWeight', 'bold');
else 
    y_txt_val= y_txt_val+2*y_txt_gap;
end

% text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('$r=%.2f$', mdl.Coefficients.Estimate(2)), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex', 'FontWeight', 'bold');
if mdl.Coefficients.pValue(2)>.05
    text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('$p=%.2f$', mdl.Coefficients.pValue(2)), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex', 'FontWeight', 'bold');
elseif mdl.Coefficients.pValue(2)>pValThresh
    text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('$p=%.3f$', mdl.Coefficients.pValue(2)), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex', 'FontWeight', 'bold');
else
    text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('$p<%.3f$', pValThresh), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex', 'FontWeight', 'bold');
end
text(x_txt_val,y_txt_val-3*y_txt_gap,sprintf('$R^2=%.2f$', mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex', 'FontWeight', 'bold');