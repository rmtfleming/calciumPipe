%**************************************************************************
%     This file is part of NeuronCyto Version 1.0.
% 
%     NeuronCyto is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     NeuronCyto is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with NeuronCyto.  If not, see <http://www.gnu.org/licenses/>.
%
%**************************************************************************
function dev=function_stdEst2D(z,method)
% Estimate noise standard deviation (AWGN model)
%
% dev = function_stdEst2D(z,method)
%
%
% OUTPUT
% ------
% dev    :  estimated noise standard deviation
%
% INPUTS
% ------
% z      :  noisy observation
% method :  method to attenuate signal (optional input)
%             0  standard shifted differences
%             1  cascaded horizontal-vertical shifted differences
%             2  wavelet domain estimation (DEFAULT)
%             3  wavelet domain estimation with boundary removal
%             4  Immerkaer's method
%
%             7  Immerkaer's method with Daubechies-based Laplacian
%             8  Blockwise Immerkaer's method with Daubechies-based Laplacian
%
%
%  methods 0-3 are based on the Median of Absolute Deviation (MAD)
%  technique, whereas method 4 is based on Laplacian filtering
%
% Alessandro Foi - Tampere University of Technology - 2005
% -----------------------------------------------------------------------

if nargin==1
    method=2; %% by default use Daubechies wavelets
end
if min(abs(method-[0 1 2 3 4 7 8]))>0
    disp(' ');disp('   !!!!!   Second argument must be either 0, 1, 2, 3, or 4.  (see help)');disp(' ');
end

if method==0 %%%% standard shifted differences
    z1=circshift(z,[1 0]);
    dz=abs(z(2:end,:)-z1(2:end,:));
    dev=median(dz(:))/(0.6745*sqrt(2));
    z2=circshift(z,[0 1]);
    dz=abs(z(:,2:end)-z1(:,2:end));
    dev=0.5*dev+0.5*median(dz(:))/(0.6745*sqrt(2));
end
if method==1   %%%  cascaded horizontal-vertical shifted differences
    z1=circshift(z,[1 0]);
    dz=z(2:end,:)-z1(2:end,:);
    dev=median(dz(:))/(0.6745*sqrt(2));
    z1=circshift(dz,[0 1]);
    dz=abs(dz(:,2:end)-z1(:,2:end));
    dev=median(dz(:))/(0.6745*2);
end

if method==2  %%% wavelet domain estimation (Daubechies length 6)
    daub6kern=[0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];
    wav_det=filter2(daub6kern,z);
    wav_det=filter2(daub6kern',wav_det);
    dev=median(abs(wav_det(:)))/.6745;
end

if method==3  %%% wavelet domain estimation (Daubechies length 6)
    daub6kern=[0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];
    wav_det=conv2(z,daub6kern,'valid');
    wav_det=conv2(wav_det,daub6kern','valid');
    dev=median(abs(wav_det(:)))/.6745;
end

if method==4    %%% Immerkaer
    LAPL=[1 -2 1;-2 4 -2;1 -2 1];
    LAPL=LAPL*sqrt(pi/2/sum(LAPL(:).^2));
    YY=conv2(z,LAPL,'valid');
    dev=mean(abs(YY(:)));
end

if method==7    %%% Immerkaer-Daubechies
    %LAPL=[1 -2 1;-2 4 -2;1 -2 1];
    daub6kern=[0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];
    LAPL=conv2(daub6kern,daub6kern);
    LAPL=conv2(LAPL,daub6kern');
    LAPL=conv2(LAPL,daub6kern');
    %    LAPL=LAPL*sqrt(pi/2/sum(LAPL(:).^2));
    %    YY=conv2(z,LAPL,'valid');
    YY=conv2(z,daub6kern*sqrt(pi/2/sum(LAPL(:).^2)),'valid');
    YY=conv2(YY,daub6kern,'valid');
    YY=conv2(YY,daub6kern','valid');
    YY=conv2(YY,daub6kern','valid');
    dev=mean(abs(YY(:)));
end


if method==8    %%% blockwise Immerkaer-Daubechies
    %LAPL=[1 -2 1;-2 4 -2;1 -2 1];
    daub6kern=[0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];
    LAPL=conv2(daub6kern,daub6kern);
    LAPL=conv2(LAPL,daub6kern');
    LAPL=conv2(LAPL,daub6kern');
    YY=conv2(z,daub6kern*sqrt(pi/2/sum(LAPL(:).^2)),'valid');
    YY=conv2(YY,daub6kern,'valid');
    YY=conv2(YY,daub6kern','valid');
    YY=conv2(YY,daub6kern','valid');
    LL=8;
    NB1=floor(size(YY,1)/LL);
    NB2=floor(size(YY,2)/LL);
    YY=YY(1:NB1*LL,1:NB2*LL);
    dev=median(mean(abs(YY(repmat(reshape(repmat([1:LL],[LL 1])'+repmat([0:LL*NB1:LL*LL*NB1-1]',[1 LL])',[LL*LL 1]) ,[1 NB1*NB2])+repmat(reshape(repmat([0:LL:NB1*LL-1],[NB2 1])+repmat(NB1*LL*[0:LL:NB2*LL-1]',[1 NB1]),[1 NB1*NB2]),[LL*LL 1]))),1));

end

