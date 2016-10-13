% Loop debugger
clear all
close all
clc

N = 100;
t = 3;

asciiProgressBar init
for count = 1:N
    asciiProgressBar( round(100.*count./N) );
%     pause( t./N );
end