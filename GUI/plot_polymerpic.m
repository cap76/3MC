function plot_polymerpic(net)

global Width

net1 = net{1};

PolymerPic(net1.Parameters,net1.X,net1.Y,net1.Z,Width);