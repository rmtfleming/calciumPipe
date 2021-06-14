function [Ucells_GT,cnt] = fuseNearbyExcitations_v3(Ucells_GT,time_fuse,t_rise,frameRate,clipping,type)

U = zeros(size(Ucells_GT));
cnt = 0;
time_fuse_frames = time_fuse*frameRate;

for i = 1:size(Ucells_GT,2)
    [imLabel,numLabels] = bwlabel(Ucells_GT(:,i));
    for j = 1:numLabels
        imLbk = (imLabel==j);
        activation = find(imLbk==1,1,'first');
        U(activation,i) = length(find(imLbk>0));
    end
    act = find(U(:,i)>0);
    j = 1;
    while j <= numLabels-1
        dist = act(j+1) - act(j);
        dur = U(act(j),i);
        separation = dist-dur;
        if separation < time_fuse_frames
            cnt = cnt+1;
            excitation_duration = dist+U(act(j+1),i); %in frames
            rand_t_rise = randi([round(t_rise(1)) round(t_rise(2))],1);
            length_rise = round((rand_t_rise/1000)*frameRate); %in frames
            decay_duration = excitation_duration - length_rise;
            decay_seconds = decay_duration/frameRate;
            newTau = decay_seconds/3;
            fused_act_pattern = generateActivationPattern(type,newTau,rand_t_rise,frameRate,clipping);
            len = length(fused_act_pattern);
            Ucells_GT((act(j):act(j)+len-1),i) = fused_act_pattern;
            if len<excitation_duration
                diference = excitation_duration - len;
                Ucells_GT((act(j)+len:act(j)+len+diference-1),i) = 0;
            end
            U(act(j),i) = len;
            U(act(j+1),i) = 0;
            act(j+1) = [];
            numLabels = numLabels-1;
            j = j-1;
        end
        j = j+1;
    end
end