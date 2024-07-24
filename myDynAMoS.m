function DynAMoS = myDynAMoS

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function fnorm = preprocess(signals,fs)
        [b,a] = butter(4,1.5/(fs/2));
        if size(signals,2)>size(signals,1)
            f_sig = filtfilt(b,a,signals');
            fnorm = vecnorm(f_sig');
        else
            f_sig = filtfilt(b,a,signals);
            fnorm = vecnorm(f_sig');
        end
    end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function [int] = intervals(act)
        % Author: Gregorio Dotti (gregorio.dotti@polito.it)
        % Last Update: 15/06/2021
        % function intervals
        %   This function takes as input the binary string which represents the
        %   muscle activation during the gait cycle and check whether it is binary,
        %   if it is normalized correctly and then extract the ON and OFF instants
        %   and the number of activations identified in a gait cycle.
        %   Input variable: act: 1x1000 binary string.
        %   Output variables: int: vector [ON1 OFF1 .... ONn OFFn], nact: number
        %   (n) equal to the number of activations.

        % check if the string is binary
        if ~isempty(find(act ~= 0 & act ~= 1))
            error('The cycles must be binary string in representation of muscle activation over the cycle')
        end

        % initialization of the output variable
        int = [];
        % identification of the samples where the muscle is active
        active = find(act);
        flag = 0;
        % identification of always off cycles
        if isempty(active)
            flag = 1;
        else
            gap = diff(active);
            ind = find(gap>1);
            % identification of the 1-madality cycles
            if isempty(ind)
                int = [active(1) active(end)];
            else
                % extraction of the activation intervals
                for i = 1:length(ind)+1
                    switch i
                        case 1
                            int = [int;active(1) active(ind(i))];
                        otherwise
                            p1 = active(ind(i-1)+1);
                            if i > length(ind)
                                p2 = active(end);
                            else
                                p2 = active(ind(i));
                            end
                            int  = [int; p1 p2];
                    end
                end
            end
        end
    end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function [ints] = ints_correction(ints,signal)
        for rep = [1:5]
            lg_ints = ints(:,2)-ints(:,1);

            flags = zeros(length(lg_ints),1);
            v_low = 0.8;
            v_up = 1.4;
            th_low = v_low*median(lg_ints);
            th_high = v_up*median(lg_ints);

            flags(lg_ints<th_low) = 1;
            flags(lg_ints>th_high) = 2;
            idxs_bef = length(find(flags));
            zeroid = 0;
            while ~isempty(find(flags))

                id = find(flags);
                id = id(1);
                if id == 1 && flags(id) == 1
                    if (ints(id+1,2)-ints(id,1))<th_high
                        ints(id,2) = ints(id+1,2);
                        ints(id+1,:) = [];
                        flags(id+1) = [];
                    else
                        flags(id) = 0;
                        zeroid = zeroid+1;
                    end

                elseif id == length(flags) && flags(id) == 1
                    if (ints(id,2)-ints(id-1,1))<th_high
                        ints(id-1,2) = ints(id,2);
                        ints(id,:) = [];
                        flags(id) = [];
                    else
                        flags(id) = 0;
                        zeroid = zeroid+1;
                    end


                else
                    switch flags(id)
                        case 1
                            if (ints(id+1,1)-ints(id,2))<(ints(id,1)-ints(id-1,2))
                                if (ints(id+1,2)-ints(id,1))<th_high
                                    ints(id,2) = ints(id+1,2);
                                    ints(id+1,:) = [];
                                    flags(id+1) = [];
                                else
                                    flags(id) = 0;
                                    zeroid = zeroid+1;
                                end
                            else
                                if (ints(id,2)-ints(id-1,1))<th_high
                                    ints(id-1,2) = ints(id,2);
                                    ints(id,:) = [];
                                    flags(id) = [];
                                else
                                    flags(id) = 0;
                                    zeroid = zeroid+1;
                                end

                            end
                        case 2
                            snippet = signal(ints(id,1):ints(id,2));
                            [~,pos] = islocalmin(snippet);
                            pos = find(pos);
                            val = snippet(pos);
                            [val,idx] = sort(val,'ascend');
                            pos = pos(idx);
                            flag = 0;
                            lmin = 1;
                            while flag ==0 && lmin<=length(val)
                                if pos(lmin) >= th_low && length(snippet)-pos(lmin)>=th_low
                                    if pos(lmin) <= th_high && length(snippet)-pos(lmin)<=th_high
                                        flag = 1;
                                        ints = cat(1,ints,[ints(id,1) ints(id,1)+pos(lmin)-1; ints(id,1)+pos(lmin)+1 ints(id,2)]);
                                        ints(id,:) =[];
                                        ints = sort(ints);
                                    else
                                        lmin = lmin+1;
                                    end
                                else
                                    lmin = lmin+1;
                                end
                            end
                            if lmin >= length(val)
                                flags(id) = 0;
                                zeroid = zeroid+1;
                            end
                            % flags(id) = 0;
                            % zeroid = zeroid+1;
                    end
                end
                lg_ints = ints(:,2)-ints(:,1);

                th_low = v_low*median(lg_ints);
                th_high = v_up*median(lg_ints);
                flags = zeros(length(lg_ints),1);
                flags(lg_ints<th_low) = 1;
                flags(lg_ints>th_high) = 2;
                idxs = find(flags);
                if (idxs_bef-length(idxs))~=0
                    idxs_bef = length(idxs);
                    zeroid = 0;
                else
                    if zeroid>0
                        flags(idxs(1:zeroid)) = 0;
                    end
                end

            end
        end
    end
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function onoff = dynamos(fnorm)

        th = max(fnorm)*0.11;
        binonoff = zeros(size(fnorm));
        binonoff(fnorm>th) = 1;
        ints = intervals(stasp);
        onoff = ints_correction(ints,fnorm);
        
    end
end