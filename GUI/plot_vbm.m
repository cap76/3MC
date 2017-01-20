function plot_vbm(net)

[D] = LDM(net);
imagesc(D);
colorbar;
