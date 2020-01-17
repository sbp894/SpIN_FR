function cMap= get_colormap(colorStart, colorEnd, nColors)
cMap(:,1)= linspace(colorStart(1), colorEnd(1), nColors);
cMap(:,2)= linspace(colorStart(2), colorEnd(2), nColors);
cMap(:,3)= linspace(colorStart(3), colorEnd(3), nColors);
end
