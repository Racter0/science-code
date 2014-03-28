function output_DQE= FindDQE(mtfmode, gain_fac, inputfile, use_noise_file, lsf_filter_level, nps_squash, no_align, divide_by_2, steep, noise_file)
%NON-INTERACTIVE VERSION
%mode 1: Eddge Spread Function
%mode 2: Envelope Function of FT

close all hidden

%numh=number of crosshatches;
%overs=oversampling factor for MTF determination
%lenh=length of crosshatches in pixels;
%samp=sampling rate of crosshatch in 1/pixel;
%nump=number of points on the crosshatch
%\
% /
%  \
numh=40;
overs=5;
%lenh needs to be uneven.
%55 worked. 45 seemed better
lenh=30*overs;
%samp=1;
samp=5;
nump=lenh/samp;
%optional turn alignment off switch
%no_align=0;

%ask for input mrc file, pixel size in um and approximate vertex locations
%if use_noise_file == 'y'
%    noise_file=input('Enter noise file name \n', 's');
%end
%also ask for vertices of a rectangular noise region.
%camera_type = input('Enter camera type: fal, u, o, k, s, f416-200, f416-120, de12 \n', 's');
%inputfile = input('Enter file name \n', 's');
    %MRC_write = input('Output MRCs? y/n \n', 's');
    %Find_DQE0 = input('Find DQE(0) independently? (takes time) y/n \n', 's');
    %filter_level = input('How much smoothing for the ESF? fit/high/med/low/none \n', 's');
filter_level = 'fit';
%lsf_filter_level =input('How much smoothing for the LSF? high/med/low/none \n', 's');
%nps_squash = input(['Sometimes NPS has low frequency artifacts. How much of the low'...
%    'end of the NPS should be flattened? \n. Give a decimal value from 0-0.1' ...
%    'where 0.1 means the first 10% of the NPS should be flattened \n']);
%divide_by_2 =input('Saved in SerialEM using 16-bit divide by 2 option? y/n \n', 's');
[map s,mi,ma,av]=ReadMRC(inputfile, 1, 1, 0);
fig=figure;
if divide_by_2 == 'y'
    disp('multiplying by 2')
    map=map*2;
end
%make square
size1=size(map, 1);
size2=size(map, 2);
if size1 > size2
    map = map(1:size2, 1:size2);
elseif size2 > size1
    map = map(1:size1, 1:size1);
end
%work with transpose
mapt=single(map)';
% switch camera_type
%     case 'u'
%         disp('Using the Ultrascan 4000')
%         gain_fac=22;
%     case 'o'
%         disp('Using the Orius')
%         gain_fac=12;
%     case 'fal'
%         disp('Using the Falcon')
%         gain_fac=189;
%     case 'fal2'
%         disp('From Henderson Falcon')
%         gain_fac=137;
%     case 'k'
%         disp('Using the K2 summit')
%         gain_fac=1;
%         %gain_fac=38.44;
%     case 's'
%         disp('Using the SC Orius')
%         gain_fac=3.2;
%     case 'f416-120'
%         disp('Using the F416')
%         gain_fac=59;
%     case 'f416-200'
%         disp('Using the F416')
%         gain_fac=36;
%     case 'de12'
%         disp('Using the de12')
%         gain_fac=30;
%     otherwise
%         disp('choose f for Falcon, o for Orius, u for Ultrascan')
%         exit
% end
mapstd=std(mapt(:))
mapmean=mean(mapt(:))
image(mapt, 'CDataMapping', 'scaled')
caxis([mapmean-mapstd, mapmean+mapstd])
%zoom on
%zoom(2)
%zoom off
colormap(gray)
set(gca, 'YDir', 'normal')
hold on

switch mtfmode
    case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Find MTF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dcm_obj = datacursormode(fig);
        set(dcm_obj,'DisplayStyle','window',...
           'SnapToDataVertex','on','Enable','on')
        disp('Select four vertices of desired area, starting with top left corner and proceeding clockwise.')
        disp('Choose an approximately rectangular area, containing a single edge, preferably half edge and half background.')
        disp('Select top left corner, then press any key')
        pause
       
        %        w0 = waitforbuttonpress;
%        w1 = waitforbuttonpress;
%        if w1 == 0
           c_info = getCursorInfo(dcm_obj);
           disp('c_info=')
           disp(c_info)
           vertex1= c_info.Position;
           scatter(vertex1(1), vertex1(2), 'r')
 %       end
        disp('Select top right corner, then press any key')
           pause
      %  w2 = waitforbuttonpress;
      %  if w2 == 0
           c_info = getCursorInfo(dcm_obj);
           disp('c_info=')
           disp(c_info)
           vertex2= c_info.Position;
           scatter(vertex2(1), vertex2(2), 'r')
      %  end
        disp('Select bottom right corner, then press any key')
       % w3 = waitforbuttonpress;
      %  if w3 == 0
           pause
           c_info = getCursorInfo(dcm_obj);
           disp('c_info=')
           disp(c_info)
           vertex3=c_info.Position;
           scatter(vertex3(1), vertex3(2), 'r')
      %  end
        disp('Select bottom left corner, then press any key')
      %  w4 = waitforbuttonpress;
      %  if w4 == 0
           pause
           c_info = getCursorInfo(dcm_obj);
           disp('c_info=')
           disp(c_info)
           vertex4=c_info.Position;
           scatter(vertex4(1), vertex4(2), 'r')
      %  end
% 	vertex1=[950, 2400];
% 	vertex2=[1200, 2400];
% 	vertex3=[1200, 2000];
% 	vertex4=[950, 2000];
        %Choose smallest possible area defined by the four corners
        chosen_rect=map(max(vertex1(1), vertex4(1)):min(vertex2(1), vertex3(1)),...
            max(vertex3(2), vertex4(2)):min(vertex1(2), vertex2(2))); 
        sc=max(vertex1(1), vertex4(1));
        uc=min(vertex2(1), vertex3(1));
        sr=max(vertex3(2), vertex4(2));
        ur=min(vertex1(2), vertex2(2));
      
        %Edge detection from rectangle:
        disp('chosen_rect mean')
        crmean=mean(chosen_rect(:))
        crstd=std(single(chosen_rect(:)))
        max(chosen_rect(:))
        min(chosen_rect(:))
  %       if   strcmp(class(chosen_rect), 'single')
  %           chosen_rect=uint16(chosen_rect);
  %           threshlevel=0.5
  %           threshlevel=graythresh(chosen_rect)
  %       else
  %           chosen_rect=uint16(chosen_rect);
  %           threshlevel=graythresh(chosen_rect)+graythresh(chosen_rect)
  %       end
%         bin_rect=im2bw(chosen_rect, threshlevel);
        bin_rect=zeros(size(chosen_rect));
        map_clean=zeros(size(map));
%	crmean=threshlevel;
        %for i=1:size(chosen_rect, 1)*size(chosen_rect, 2)
        for i=1:size(map, 1)*size(map, 2) 
            if map(i)>crmean
                map_clean(i)=1;
            else
                map_clean(i)=0;
            end
        end

         %Go through and delete single specks of noise (identifiable as 010 or 101)
         three_mat=zeros(1, 3);

        for i=1:size(map, 2)
            for j=3:size(map, 1)
                three_mat(1)=three_mat(2);
                three_mat(2)=three_mat(3);
                three_mat(3)=map_clean(j, i);
                if (three_mat == [0, 1, 0])
                    map_clean(j-1, i)=0;
                end
                if (three_mat == [1, 0, 1])
                    map_clean(j-1, i)=1;
                end
            end
        end

        three_mat=[0,0,0];
        for j=1:size(map, 1)
            for i=3:size(map, 2)
                three_mat(1)=three_mat(2);
                three_mat(2)=three_mat(3);
                three_mat(3)=map_clean(j, i);
                if (three_mat == [0, 1, 0])
                    map_clean(j, i-1)=0;
                end
                if (three_mat == [1, 0, 1])
                    map_clean(j, i-1)=1;
                end
            end
        end

        %we have so far been altering the entire map; now choose only the
        %selected rectange.
        submap=map_clean(max(vertex1(1), vertex4(1)):min(vertex2(1), vertex3(1)),...
            max(vertex3(2), vertex4(2)):min(vertex1(2), vertex2(2))); 
        figure;
        subplot(2, 2, 1), image(submap, 'CDataMapping', 'scaled')
        caxis([0 1])
        colormap(gray)
        set(gca, 'YDir', 'normal')
        subplot(2, 2, 2), image(chosen_rect, 'CDataMapping', 'scaled')
        caxis([mapmean-mapstd, mapmean+mapstd])
        colormap(gray)
        set(gca, 'YDir', 'normal')
        subplot(2, 2, 3), image(map_clean, 'CDataMapping', 'scaled')
        caxis([0 1])
        colormap(gray)
        set(gca, 'YDir', 'normal')
        subplot(2, 2, 4), image(map, 'CDataMapping', 'scaled')
        caxis([mapmean-mapstd, mapmean+mapstd])
        colormap(gray)
        set(gca, 'YDir', 'normal')
        diffmat1=abs(diff(submap, 1, 1));
        diffmat2=abs(diff(submap, 1, 2));
        diffmat=zeros(size(submap, 1)-1, size(submap, 2)-1);
        for i=1:(size(submap, 1)-1)
            for j=1:(size(submap, 2)-1)
                if diffmat1(i, j)==1 || diffmat2(i, j)==1
                    diffmat(i, j)=1;
                end
            end
        end
        %working with transpose due to display conventions
        %at the moment it is not a transpose. Let's see!
        diffmatt=diffmat;
        figure;
        imshow(diffmatt)
        title('diffmatt')
        %denoise by finding and destroying small clusters
        ncounts=0;
        for i=1:size(diffmatt, 1)
            for j=1:size(diffmatt, 2)
                if diffmatt(i, j)==1
                    for k=-7:7
                        for l=-7:7
                            if i+k > 0 && i+k < size(diffmatt, 1) && j+l >0 && j+l < size(diffmatt, 2)
                                if diffmatt(i+k, j+l) == 1
                                    ncounts=ncounts+1;
                                end
                            end
                        end
                    end
                    if ncounts <= 7
                        diffmatt(i, j)=0;
                    end
                    ncounts=0;
                end
            end
        end
        figure;
        imshow(diffmatt)
        title('diffmatt')
        
        
        %When the line is very vertical, the diagonal is represented by an
        %insufficient # of points. Determine roughly whether the slope is
        %more vertical or horizontal.
        %look for first point on diag in first row
        [C, x2]=max(diffmatt(1, :));
        y2=1;
        if C == 0
            disp('point 1')
            %then try in last column
            [C, y2]=max(diffmatt(:, size(diffmatt, 2)));
            x2=size(diffmatt, 2);
        end
        %Then look in last row
        [C, x1]=max(diffmatt(size(diffmatt, 1), :));
        y1=size(diffmatt, 1);
        if C == 0
            disp('point 2')
            %or first column
            [C, y1]=max(diffmatt(:, 1));
            x1=1;
        end
        disp('x2-x1')
        abs(x2-x1)
        disp('y2-y1')
        abs(y2-y1)
        
   %     if abs(x2-x1) > abs(y2-y1)
            %Then line is more horizontal
            diagonal=zeros(size(diffmatt, 1), 1);
                for i=1:size(diffmatt, 1)
                    [C, diagonal(i)]=max(diffmatt(i, :));
                end
  %      else
  %          %line is more vertical
  %          diagonal=zeros(size(diffmatt, 2), 1);
  %              for i=1:size(diffmatt, 2)
  %                  [C, diagonal(i)]=max(diffmatt(:, i));
  %              end
  %      end
        %to not count the flat parts of a diagonal
        zcount=0;
        buf1=0;
        setbuf=1;
        for i=1:length(diagonal)
            if diagonal(i) <=10
                zcount=zcount+1;
            else
                if setbuf==1
                    buf1=zcount;
                    setbuf=0;
                end
                zcount=0;
            end
        end
        diagonal_new=diagonal(1+buf1:(length(diagonal)-zcount));
        ld=length(diagonal_new);
        x=(1:1:ld)';
        p=polyfit(x, diagonal_new, 1);
        %measure slope angle
        hyp=sqrt((p(1)*(ld+buf1)+p(2)-(p(1)*(1+buf1)+p(2)))^2+(ld-1)^2);
        opp=p(1)*(ld+buf1)+p(2)-(p(1)*(1+buf1)+p(2));
        angle=asin(opp/hyp)*(360/(2*pi))
        figure;
        subplot(2, 2, 1), image(diffmatt', 'CDataMapping', 'Scaled')
        set(gca, 'YDir', 'normal')
        axis equal
        colormap(gray)
        subplot(2, 2, 2), image(chosen_rect', 'CDataMapping', 'Scaled')
        hold on
        set(gca, 'YDir', 'normal')
        axis equal
        colormap(gray)
        scatter(x+buf1, diagonal_new, 'o')
        subplot(2, 2, 3), image(chosen_rect', 'CDataMapping', 'Scaled')
        set(gca, 'YDir', 'normal')
        axis equal
        colormap(gray)
        hold on
        plot(x+buf1, (p(1)*(x)+p(2)), 'o')

        
        %nh=length of diagonal
        nh=floor(sqrt((p(1)*x(1)-p(1)*x(ld))^2+(ld-1)^2));
        disp('components of nh')
        sqrt((p(1)*x(1)-p(1)*x(ld))^2)
        sqrt((ld-1)^2)
        %make perpendicular crosshatches
        crosshatches=zeros(lenh, ld);
        crosshatch_bin=zeros(lenh, ld);
        hatch_coordsr=zeros(lenh, ld);
        hatch_coordsc=zeros(lenh, ld);
        mapt_bin=map_clean';
        %yint=yshift for a diagonal of 1 pixel
        %overs= overasmpling factor.
        theta=asin((p(1)*x(ld))/nh);
        disp('p(1)');
        p(1);
        x(ld);
        nh;
        gamma=pi/2-theta;
        yint=sin(gamma)/overs;
        xint=cos(gamma)/overs;
        for i=1:lenh
            for j= 1:ld
                crosshatches(i, j)=mapt(ceil(sr+p(1)*x(j)+p(2)+yint*(lenh/2+1-(i-1))), ceil(buf1+x(j)+sc-xint*(lenh/2+1-(i-1))));
                crosshatch_bin(i, j)=mapt_bin(ceil(sr+p(1)*x(j)+p(2)+yint*(lenh/2+1-(i-1))), ceil(buf1+x(j)+sc-xint*(lenh/2+1-(i-1))));
                hatch_coordsr(i, j)=ceil(sr+p(1)*x(j)+p(2)+yint*(lenh/2+1-(i-1)));
                hatch_coordsc(i, j)=ceil(buf1+sc+x(j)-xint*(lenh/2+1-(i-1)));
            end
        end
        subplot(2, 2, 4), image(mapt, 'CDataMapping', 'Scaled')
        hold on
        set(gca, 'YDir', 'normal')
        axis equal
        caxis([mapmean-mapstd, mapmean+mapstd])
        colormap(gray)
        scatter(hatch_coordsc(:), hatch_coordsr(:))
        figure;
        image(crosshatches, 'CDataMapping', 'Scaled')
        set(gca, 'YDir', 'normal')
        axis equal
        title('crosshatches')
        colormap(gray)
        figure;
        image(crosshatch_bin, 'CDataMapping', 'Scaled')
        set(gca, 'YDir', 'normal')
        axis equal
        title('crosshatches-bin')
        colormap(gray)
        %align the crosshatches to make up for wavy edge.
        %lags are shifting b rightwards relative to a.

        %try making crosshatches more alignable.
%        av_value=mean(crosshatches);
%        crosshatch_bin=zeros(size(crosshatches));
%         for i=1:size(crosshatches, 1)
%             for j=1:size(crosshatches, 2)
%                 if (crosshatches(i, j) >= av_value)
%                     crosshatch_bin(i, j)=1;
%                 end
%             end
%         end
        
        %%Added in some -3 and -6 in order to get rid of the edges which are still
        %%noisy

        %%%%correlate derivative

%         deriv_hatches1=abs(diff(crosshatch_bin, 1, 1));
%         deriv_hatches2=abs(diff(crosshatch_bin, 1, 2));
%         deriv_hatches=zeros(size(deriv_hatches1, 1), size(deriv_hatches2, 2));
%         for i=1:size(deriv_hatches, 1)
%             for j=1:size(deriv_hatches, 2)
%                 if deriv_hatches1(i, j) == 1 || deriv_hatches2(i, j)==1
%                     deriv_hatches(i, j)=1;
%                 end
%             end
%         end

        %Only take derivative along one dimension because otherwise there
        %is a false extra pixel of thickness to bumps
        deriv_hatches=abs(diff(crosshatch_bin, 1, 1));
        
        figure;
        image(deriv_hatches', 'CDataMapping', 'Scaled')
        set(gca, 'YDir', 'normal')
        axis equal
        title('deriv_hatches')
        colormap(gray)
        
       
        
%        if abs(angle) <= 45
            lag_mat=zeros(size(deriv_hatches, 2), 1);
            size(deriv_hatches)
            middlehatch=ceil(size(deriv_hatches, 2)/2);
            %deriv_hatches(:, middlehatch)
        for i=1:size(deriv_hatches, 2)
            [c1, lags]=xcorr(deriv_hatches(:, middlehatch), deriv_hatches(:,i));
            [C, I]=max(c1);
            lag_mat(i)=lags(I);
        end
        lag_mat
%         elseif abs(angle) > 45
%             lag_mat=zeros(size(deriv_hatches, 1), 1);
%             size(deriv_hatches)
%             middlehatch=ceil(size(deriv_hatches, 2)/1);
%              for i=1:size(deriv_hatches, 1)
%                 [c1, lags]=xcorr(deriv_hatches(middlehatch, :), deriv_hatches(i, :));
%                 [C, I]=max(c1);
%                 lag_mat(i)=lags(I);
%              end
%         end

        %Note: the deriv_hatches are 1 smaller than the crosshatches
        maxdiff=(max(lag_mat)-min(lag_mat));
        %crosshatch_ali=zeros(nump+maxdiff, ld-6);
        crosshatch_ali=zeros(size(crosshatches)-[1, 1]);
 %       if abs(angle) <=45
            for j=1:ld-6
                for i=1:lenh
                    crosshatch_ali(i-min(lag_mat)+lag_mat(j), j)=crosshatches(i,j+3);
                end
            end
            averagech_long=sum(crosshatch_ali, 2)/(ld-6);
            averagech=averagech_long(1+maxdiff:length(averagech_long)-maxdiff);
         
%         elseif abs(angle) > 45
%             for j=1:size(crosshatches, 1)-1 %lenh
%                 for i=1:size(crosshatches, 2)-1 %ld-6
%                     crosshatch_ali(j, i-min(lag_mat)+lag_mat(j))=crosshatches(j, i);
%                 end
%             end
%             averagech_long=sum(crosshatch_ali, 1)/(size(crosshatches, 1)-1);
%             disp('size(averagrech_long')
%             %trim the  alignment bits (filled with zeros) off
%             averagech=averagech_long(1, 1+maxdiff:length(averagech_long)-maxdiff);
%             size(averagech_long)
%         end
        esf=averagech;
        if size(esf, 1) < size(esf, 2)
            esf=esf';
        end
        disp('size esf')
        size(esf)
        if no_align
            esf=mean(crosshatches, 2);
        end
        
        if esf(1) > esf(length(esf))
            esf=flipud(esf);
        end
        x=(1: 1: length(esf))'/overs;
        
         figure;
        image(crosshatch_ali, 'CDataMapping', 'scaled')
        colormap(gray)
        axis equal
        title('aligned edge')
        set(gca, 'YDir', 'normal')      
        
        
        %fit double gaussian to ESF as in McMullan 2009
%         xcent=133/5;
%         dblgaussian=@(lambda, xarg) (lambda(2)*lambda(1)/2)*(((((xarg-xcent)+0.5)/lambda(1)).* ...
%             erfc((-(xarg-xcent)-0.5)/lambda(1))-1/sqrt(pi)*exp(-(((xarg-xcent)+0.5)/lambda(1)).^2) ...
%             - (((xarg-xcent)-0.5)/lambda(1)).*erfc((-(xarg-xcent)+0.5)/lambda(1))+1/sqrt(pi) ...
%             *exp(-(((xarg-xcent)-0.5)/lambda(1)).^2)));
%         fitparams0=[1, mean(mean(map(1:500, 1:500)))*0.5]
%         %figure;
%         %plot(x, dblgaussian(fitparams0, x))
%         fitparams=nlinfit(x, esf, dblgaussian, fitparams0)
%         %figure;

        
        switch filter_level
            case 'fit'
                %curve fit instead of smoothing
                %figure out center for knot placement
                [C, cent]=max(diff(esf));
                %esf_slm = slmengine(x, esf, 'knots', [min(x), 5, ...
                %    (cent/overs)-6, (cent/overs)-4, (cent/overs)-3, (cent/overs)-2:0.5:(cent/overs)+3,...
                %    (cent/overs)+4, 23, max(x)],'plot', 'on');
                esf_slm = slmengine(x, esf, 'knots', [min(x), 5, ...
                    (cent/overs)-4, (cent/overs)-3:0.5:(cent/overs)+3,...
                    (cent/overs)+4, 23, max(x)],'plot', 'on', 'increasing', 'on', 'concaveup', [cent/overs-4, cent/overs],...
                    'concavedown', [cent/overs, cent/overs+4])
                if steep == 'y'
                    esf_slm = slmengine(x, esf, 'knots', [min(x), 5, ...
                    (cent/overs)-4, (cent/overs)-3:0.5:(cent/overs)+3,...
                    (cent/overs)+4, 23, max(x)],'plot', 'on')
                end
                esf_smooth=slmeval(x, esf_slm, 0);
            case 'high'
                 %smooth ESF as in Samei 1998
                window=17;
                pweights=zeros(window, 1);
                esf_smooth=esf;
                for i=-(window-1)/2:(window-1)/2
                    pweights(i+1+(window-1)/2)=exp(-(4*i/(window-1))^2);
                end
                for i=1:length(esf)-(window+1)
                    [p, s, mu]=wpolyfit(x(i:i+window-1), esf_smooth(i:i+window-1), 2, pweights);
                    xi=(window-1)/2+i;
                    esf_smooth(xi)=polyval(p, x(xi), s, mu);
                end
            case 'med'
                window=17;
                pweights=zeros(window, 1);
                esf_smooth=esf;
                for i=-(window-1)/2:(window-1)/2
                    pweights(i+1+(window-1)/2)=exp(-(4*i/(window-1))^2);
                end
                for i=1:length(esf)-(window+1)
                    [p, s, mu]=wpolyfit(x(i:i+window-1), esf_smooth(i:i+window-1), 4, pweights);
                    xi=(window-1)/2+i;
                    esf_smooth(xi)=polyval(p, x(xi), s, mu);
                end
            case 'low'
                esf_smooth=esf;
            case 'none'
                esf_smooth=esf;
            otherwise %same as 'med'
                window=17;
                pweights=zeros(window, 1);
                esf_smooth=esf;
                for i=-(window-1)/2:(window-1)/2
                    pweights(i+1+(window-1)/2)=exp(-(4*i/(window-1))^2);
                end
                for i=1:length(esf)-(window+1)
                    [p, s, mu]=wpolyfit(x(i:i+window-1), esf_smooth(i:i+window-1), 2, pweights);
                    xi=(window-1)/2+i;
                    esf_smooth(xi)=polyval(p, x(xi), s, mu);
                end
        end
        

        %Differentiate ESF to get LSF; smooth with Hamming filter (Samei
        %1998)
        lsf=diff(esf_smooth)/(1/sqrt(2)*overs);
        switch lsf_filter_level
            case 'fit'
                %no additional smoothing should be required
                lsf_hammed=lsf;
            case 'high'
                ham_window=zeros(length(lsf), 1);
                ham_window_curve=hamming(ceil(length(lsf)*1/3));
                ham_window(1+ceil(length(ham_window)*1/3):length(ham_window_curve)+ceil(length(ham_window)*1/3))=ham_window_curve;
                lsf_ft=fftshift(fft(lsf));
                lsf_hammed1=ifft(fftshift(lsf_ft.*ham_window));
                lsf_hammed=sqrt(real(lsf_hammed1).^2+imag(lsf_hammed1).^2);
            case 'med'
                ham_window=zeros(length(lsf), 1);
                ham_window_curve=hamming(ceil(length(lsf)*3/5));
                ham_window(1+ceil(length(ham_window)*1/5):length(ham_window_curve)+ceil(length(ham_window)*1/5))=ham_window_curve;
                lsf_ft=fftshift(fft(lsf));
                lsf_hammed1=ifft(fftshift(lsf_ft.*ham_window));
                lsf_hammed=sqrt(real(lsf_hammed1).^2+imag(lsf_hammed1).^2);
            case 'low'
                ham_window=hamming(length(lsf));
                lsf_ft=fftshift(fft(lsf));
                lsf_hammed1=ifft(fftshift(lsf_ft.*ham_window));
                lsf_hammed=sqrt(real(lsf_hammed1).^2+imag(lsf_hammed1).^2);
            case 'none'
                lsf_hammed=lsf;
            otherwise % same as med
                ham_window=zeros(length(lsf), 1);
                ham_window_curve=hamming(ceil(length(lsf)*3/5));
                ham_window(1+ceil(length(ham_window)*1/5):length(ham_window_curve)+ceil(length(ham_window)*1/5))=ham_window_curve;
                lsf_ft=fftshift(fft(lsf));
                lsf_hammed1=ifft(fftshift(lsf_ft.*ham_window));
                lsf_hammed=sqrt(real(lsf_hammed1).^2+imag(lsf_hammed1).^2);
        end
        if lsf_hammed(ceil(length(lsf_hammed)/2+1)) < 0
            disp('flipping lsf_hammed')
            lsf_hammed=-lsf_hammed;
        end

        
        lsf_pad=padarray(lsf, 4*length(lsf));
        lsf_hammed_pad=padarray(lsf_hammed, 4*length(lsf_hammed));
        %spacing for abcissa: 1 / (N dt)
        % where N is # of samples and dt is sampling period
        % here dt is psize*samp, and N is nump
        % abcissa at kth point is f_k=k df where df= 1 / (N dt)
        % http://www.av8n.com/physics/fourier-refined.htm
        mtf=sqrt(real(fftshift(fft(lsf_pad))).^2+imag(fftshift(fft(lsf_pad))).^2);
        mtf_from_hammed=sqrt(real(fftshift(fft(lsf_hammed_pad))).^2+imag(fftshift(fft(lsf_hammed_pad))).^2);
        disp('length of mtf_from_hammed')
        length(mtf_from_hammed)
        [mtf_fh_max, mtf_fh_max_i]=max(mtf_from_hammed);
        mtf_from_hammed=mtf_from_hammed/mtf_fh_max;
        mtf_fh=mtf_from_hammed(mtf_fh_max_i:length(mtf_from_hammed));
        
        mtf_fh_samp=mtf_fh(1:ceil(length(mtf_fh)/overs));
        disp('length of mtf_fh_samp')
        disp(length(mtf_fh_samp))
        x_fh=(0:1:(length(mtf_fh_samp)-1))'/(length(mtf_fh_samp)-1);
        [mtf_max, mtf_max_i]=max(mtf);
        mtf=mtf/mtf_max;
        mtf_0_1=mtf(mtf_max_i:length(mtf));
        mtf_0_1_filt=mtf_0_1;
        for i=ceil(length(mtf)*0.1):length(mtf_0_1)-3
            mtf_0_1_filt(i)=(mtf_0_1(i)+mtf_0_1(i+1)+mtf_0_1(i+2))/3;
        end
        mtf_smooth2=mtf_0_1_filt;
        mtf_smooth=medfilt1(mtf_0_1_filt, 6);
        mtf_smooth2(ceil(length(mtf_smooth)/3):length(mtf_smooth))=mtf_smooth(ceil(length(mtf_smooth)/3):length(mtf_smooth));
        xfreq_h=(0:1:length(mtf_0_1_filt)-1)'/length(mtf_0_1_filt);
        mtf_samp=mtf_smooth2(1:ceil(length(mtf_smooth2)/overs));
        if length(mtf_samp) > length(xfreq_h)
            mtf_samp=mtf_samp(1:length(mtf_samp)-1);
        end
       
        
        mtf_slm = slmengine(x_fh, mtf_fh_samp,'plot', 'on', 'xy', [0, 1], 'knots', 8);
        mtf=slmeval(x_fh, mtf_slm, 0);
        xfreq=(0:1:(length(mtf)-1))'/(length(mtf)-1);
        
        disp('size of mtf_fh_samp')
        size(mtf_fh_samp)
        disp('size of mtf')
        size(mtf)
        disp('size of mtf_samp')
        size(mtf_samp)
        
        figure;
        subplot(2,2,1), plot(x, esf, 'b')
        hold on
        plot(x, esf_smooth, 'g')
        hold off
        title('edge spread fn')
        ylabel('counts')
        xlabel('microns')
        
        subplot(2,2,2), plot(x(1:length(x)-1), lsf)
        hold on
        plot(x(1:length(x)-1), lsf_hammed, 'g')
        hold off
        title('line spread fn')
        subplot(2, 2, 3), plot(xfreq_h, mtf_0_1, 'r');
        plot( (0:1:(length(mtf_samp)-1))'/(length(mtf_samp)-1), mtf_samp, 'b') 
        title('mtf, mtf=r, mtf_fh_samp=g, mtf_samp=b')
        xlabel('fraction of nyquist')
        axis([0 1 -Inf Inf])
        hold on
        disp('size of xfreq')
        size(xfreq)
        plot(xfreq, mtf, 'r')
        plot(x_fh, mtf_fh_samp, 'g')
        hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Noise Power Spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

%Divide image into npatches patches of size spatches, marking the center of
%each patch with an areadot
        spatches=512; %Must be even
        npatches=floor(size(map, 1)/spatches)*floor(size(map, 2)/spatches);
        noise_areas=zeros(spatches, spatches, npatches);
        area_dots=zeros(npatches, 2);
        noise_fts=zeros(spatches, spatches, npatches);
        %padsq_length=length(fftshift(fft2(single(noise_areas(:,:, 1)), spatches*4, spatches*4)))
        %noise_fts_pad=zeros(padsq_length, padsq_length, npatches);
        nps_stack=zeros(spatches, spatches, npatches);
        %nps_stack_pad=zeros(padsq_length, padsq_length, npatches);
        noise_stats=zeros(npatches, 2);
        noise_include=ones(npatches, 1);
        index_convert=zeros(npatches, 2);
        
        %Create a cose-edge mask to go around all of the square patches.
        cew=10;
        cos_x=(0:(pi/cew):pi);
        cos_edge=1-((1+cos(cos_x))/2);
        cos_mask=ones(spatches);
        for i=1:cew
            for j=1:spatches
            if cos_mask(i, j) > cos_edge(i)
                cos_mask(i, j)=cos_edge(i);
                cos_mask(j, i)=cos_edge(i);
                cos_mask(spatches-i, j)=cos_edge(i);
                cos_mask(j, spatches-i)=cos_edge(i);
            end
            end
        end
            
           
        mask_weight=sum(sum(cos_mask.^2));
        size(map)
        for i=1:floor(size(map, 1)/spatches)
            for j=1:floor(size(map, 2)/spatches)
                index1d=(i-1)*floor(size(map, 1)/spatches)+j;
                index_convert(index1d, 1)=j;
                index_convert(index1d, 2)=i;
                noise_areas(:,:,index1d)=map((i-1)*spatches+1:(i-1)*spatches+spatches, ...
                    (j-1)*spatches+1:(j-1)*spatches+spatches, 1);
                noise_mean_tmp=mean(mean(noise_areas(:, :, index1d)));
                noise_areas(:, :, index1d)=noise_areas(:, :, index1d).*cos_mask+noise_mean_tmp*(1-cos_mask);
                area_dots(index1d, 1)=spatches*(j-1)+spatches/2;                
                area_dots(index1d, 2)=spatches*(i-1)+spatches/2;
                noise_fts(:,:,index1d)=fftshift(fft2(single(noise_areas(:,:,index1d))));
              %  noise_fts_pad(:, :, index1d)=fftshift(fft2(single(noise_areas(:,:, index1d)), spatches*4, spatches*4));         
                nps_stack(:, :, index1d)=(abs(noise_fts(:,:,index1d)).^2)/(spatches*spatches);
               % nps_stack_pad(:, :, index1d)=(abs(noise_fts_pad(:,:,index1d)).^2)/(spatches^2);
                noise_tmp=noise_areas(:, :, index1d);
                noise_stats(index1d, 1)=mean(noise_tmp(:));
                noise_stats(index1d, 2)=var(noise_tmp(:));                
            end
        end

        %subplot(2, 2, 4), image(log(nps_stack(:, :, 3)), 'CDataMapping', 'scaled')
        %colormap(gray)

%Determine which patches contain pointer by finding patches with 
%significantly lower mean values. Do statistics on only largest half of values
%to avoid pointer. Mark pointer patches with a red dot and set noise_include to 0.
        noise_mean_sort=sort(noise_stats(:, 1), 'descend');
        noise_mean_mean=mean(noise_mean_sort(1:ceil(length(noise_mean_sort)/2)));
        noise_mean_stdev=std(noise_mean_sort(1:ceil(length(noise_mean_sort)/2)));
        red_dots=0;
        for i=1:npatches
            if noise_stats(i, 1)<=noise_mean_mean-3*noise_mean_stdev
                noise_include(i)=0;
                red_dots=red_dots+1;
            end
        end
        red_dot_coords=zeros(red_dots, 2);
        j=1;
        for i=1:npatches
            if noise_include(i) == 0
                red_dot_coords(j, 1) = area_dots(i, 1);
                red_dot_coords(j, 2) = area_dots(i, 2);
                j= j+1;
            end
        end    
        figure;
        image(mapt, 'CDataMapping', 'scaled')
        colormap(gray)
        axis equal
        set(gca, 'YDir', 'normal')
        title('red dots on the pointer, blue dots on noise patches')
        hold on
        scatter(area_dots(:, 1), area_dots(:, 2))
        scatter(red_dot_coords(:, 2), red_dot_coords(:, 1), 'r')
%Average NPS; average noise value
        tmpsum_ps=zeros(size(nps_stack(:,:,1)));
        %tmpsum_ps=zeros(size(nps_stack_pad(:, :, 1)));
        tmpsum_real=zeros(size(noise_areas(:,:,1)));
        tmpsum_squared=zeros(size(noise_areas(:,:,1)));
        tmpcnts=0;
        for i=1:length(noise_include)
            if noise_include(i) == 1
                tmpsum_ps=tmpsum_ps+nps_stack(:,:,i);
%                   figure;
%                   image(log(nps_stack(:, :, i)), 'CDataMapping', 'scaled')
%                   colormap(gray)
                %tmpsum_ps=tmpsum_ps+nps_stack_pad(:,:,1);
                tmpsum_real=tmpsum_real+noise_areas(:,:,i);
                tmpsum_squared=tmpsum_squared+noise_areas(:,:,i).^2;
                tmpcnts=tmpcnts+1;
            end
        end
        nps=tmpsum_ps/tmpcnts;
        noise_squared=sum(tmpsum_squared(:))/tmpcnts/length(noise_areas)^2;
        noise_mean=mean(mean(tmpsum_real/tmpcnts))
        noise_rms=sqrt(noise_squared-noise_mean^2)
%For use in other parts of the code, make a map with all of the pointer
%squares replaced by the background value.
%         map_meanvals=double(map);
%         cew=30;
%         cos_x=(0:(pi/cew):pi);
%         cos_edge=(1+cos(cos_x))/2;
%         for i=1:length(red_dot_coords)         
%             map_meanvals(red_dot_coords(i, 2)-spatches/2+1:red_dot_coords(i, 2)+spatches/2, ...
%                 red_dot_coords(i, 1)-spatches/2+1:red_dot_coords(i, 1)+spatches/2)=noise_mean;
% %   BETTER IDEA--you should put cosedge masks centered on the dots.
%          
%         end
%        new_size=floor(length(map_meanvals)/spatches);
%        map_meanvals=map_meanvals(1:new_size*spatches, 1:new_size*spatches);
%        figure;
%        image(map_meanvals', 'CDataMapping', 'scaled')
%        title('map_meanvals')
%        colormap(gray)
       
%Radial average and plot
% 
%         tmp=fftshift(fft2(single(map(2500:3524, 2500:3524))));
%         nps=(abs(tmp).^2)/(1024*1024);
%         figure;
%         image(log(nps), 'CDataMapping', 'scaled');
%         colormap(gray);
        
        enditer=spatches^2;
        cr=spatches/2+1;
        cc=spatches/2+1;
        counts=zeros(spatches/2+1, 1);
        ps=zeros(spatches/2+1, 1);
        for i=1:enditer
            [r c]=ind2sub([spatches, spatches], i);
            radius=ceil(sqrt((r-cr)^2+(c-cc)^2));
            if radius <= length(ps)-1 && radius > 0
                counts(radius)=counts(radius)+1;
                ps(radius)=ps(radius)+nps(r, c);
            end
        end

        

     
        %Alternative way of finding the noise power spectrum: find from
        %overlapping patches in the top left 1024x1024 area
%         submap=map(2500:3524, 2500:3524);
%         area_count=0;
%         sqsize=32;
%         cew=6;
%         cos_x=(0:(pi/cew):pi);
%         cos_edge=1-((1+cos(cos_x))/2);
%         cos_mask=ones(sqsize);
%         size(cos_mask)
%         for i=1:cew
%             for j=1:sqsize
%                 if cos_mask(i, j) > cos_edge(i)
%                     cos_mask(i, j)=cos_edge(i);
%                     cos_mask(j, i)=cos_edge(i);
%                     cos_mask(sqsize-i, j)=cos_edge(i);
%                     cos_mask(j, sqsize-i)=cos_edge(i);
%                 end
%             end
%         end
%         noise_ps_cm=zeros(sqsize);
%         for i=sqsize/2+1:3:1024-sqsize/2+1
%             for j=sqsize/2+1:3:1024-sqsize/2+1
%                 area_count=area_count+1;
%                 noise_area=submap(i-sqsize/2:i+sqsize/2-1, j-sqsize/2:j+sqsize/2-1);
%                 noise_ft=fftshift(fft2(single(noise_area).*cos_mask));
%                 noise_ps=(abs(noise_ft).^2)/(sqsize*sqsize);
% %                 figure;
% %                 image(log(noise_ps), 'CDataMapping', 'scaled')
% %                 colormap(gray)
%                 noise_ps_cm=noise_ps_cm+noise_ps;
%             end
%         end
%         noise_ps_cm=noise_ps_cm/area_count;
%         figure;
%         image(log(noise_ps_cm), 'CDataMapping', 'scaled')
%         colormap(gray)
%         title('cumulative nps')
%         disp('acquired cumulative noise power spectrum')
%         %radial average and plot
%         cr=sqsize/2+1;
%         cc=sqsize/2+1;
%         counts=zeros(sqsize/2+1, 1);
%         ps_1=zeros(sqsize/2+1, 1);
%         for i=1:sqsize^2
%             [r c]=ind2sub([sqsize, sqsize], i);
%             radius=ceil(sqrt((r-cr)^2+(c-cc)^2));
%             if radius <= length(ps_1)-1 && radius > 0
%                 counts(radius)=counts(radius)+1;
%                 ps_1(radius)=ps_1(radius)+noise_ps_cm(r, c);
%             end
%         end
%         normps=bsxfun(@rdivide, ps_1, counts);

        normps=bsxfun(@rdivide, ps, counts);
        normps(length(normps))=normps(length(normps)-2);
        normps(length(normps)-1)=normps(length(normps)-2);
        
        normps_hammed=normps;
         ham_x=(0:1:length(normps_hammed)-1)/(length(normps_hammed)-1);

        flatl=ceil(length(normps_hammed)*nps_squash)
        normps_hammed(1: flatl)=mean(normps_hammed(flatl+1:flatl+ceil(flatl/3)));
      
        %Smooth the power spectrum
        normps_smooth=normps_hammed;
         for i=1:length(normps_hammed)-6
             normps_smooth(i)=(normps_hammed(i)+normps_hammed(i+1)+normps_hammed(i+2)+normps_hammed(i+3)+normps_hammed(i+4)+normps_hammed(i+5))/6;
         end
        
         figure;
         plot(ham_x, normps_smooth, 'b', ham_x, normps_hammed, 'g')
         title('smoothed=b, unsmoothed=g')
         
        nps_slm = slmengine(ham_x, normps_smooth, 'plot', 'on', 'knots', 10);
        nps_slm_eval=slmeval(ham_x, nps_slm, 0);
        normps_filt=medfilt1(nps_slm_eval, 6);
        
        %resize nps by resampling so it is the same as the mtf
        if(length(normps_filt) >= length(mtf_samp)) 
            disp('the resampling is going from')
            disp('checkpoint')
            length(normps_filt)
            length(mtf_samp)
            normps_red=resample(normps_filt, length(xfreq), length(normps_filt), 0)';
        else
            disp('the resampling is going from')
            disp('checkpoint 2')
            size(mtf_samp)
            size(normps_filt)
            x_samp=(0:1:length(normps_filt)-1)/(length(normps_filt)-1);
            x_resamp=(0:1:length(mtf_samp)-1)/(length(mtf_samp)-1);
            normps_red=interp1(x_samp, normps_filt, x_resamp)';
        end
        
        normps_red(1)=normps_red(2);
        freq_dqe=(0:1:length(normps_red)-1)/(length(normps_red)-1);
        freq_dqe_long=(0:1:length(normps_filt)-1)/length(normps_filt);
        freq_dqe_smooth=(0:1:length(normps_smooth)-1)/length(normps_smooth);
        
        figure;
        plot(freq_dqe_long, log(normps));
        hold on
        title('log(Noise Power Spectrum)')
        plot(freq_dqe_long, log(normps_filt), 'r');
        %plot(ham_x, log(normps_hammed), 'g');
        plot(freq_dqe, log(normps_red), 'or')
        plot(freq_dqe_smooth, log(normps_smooth), 'og')
        figure;
        plot(freq_dqe_smooth, normps_smooth/(normps_smooth(3)))
        title('NPS normalized to 1')
        
% The two noise images method.
    if use_noise_file == 'y'
        [noise_map, one, two, three, four]=ReadMRC(noise_file, 1, 1, 0);
        sqsize=1024;
        noise_map=noise_map';
        noise_map=single(noise_map(1:sqsize, 1:sqsize));
        class(noise_map)
        class(mapt)
        noise_diff=mapt(1:sqsize, 1:sqsize)-noise_map(1:sqsize, 1:sqsize);
        noise_ft=fftshift(fft2(single(noise_diff)));%.*cos_mask));
        noise_ps=(abs(noise_ft).^2)/(sqsize*sqsize);
         %radial average and plot
        cr=sqsize/2+1;
        cc=sqsize/2+1;
        counts=zeros(sqsize/2+1, 1);
        ps=zeros(sqsize/2+1, 1);
        for i=1:sqsize^2
            [r c]=ind2sub([sqsize, sqsize], i);
            radius=ceil(sqrt((r-cr)^2+(c-cc)^2));
            if radius <= length(ps)-1 && radius > 0
                counts(radius)=counts(radius)+1;
                ps(radius)=ps(radius)+noise_ps(r, c);
            end
        end
        %multiply by 1/sqrt(2) because of statistics
        normps_diff=bsxfun(@rdivide, ps, counts)*1/sqrt(2);
        normps_diff_x=(0:1:length(normps_diff)-1)'/(length(normps_diff)-1);
    
        %%%%%ADD MORE SMOOTHING FIRST HERE
        normps_diff_smooth=normps_diff;
        for i=1:length(normps_diff)-6
            normps_diff_smooth(i)=(normps_diff(i)+normps_diff(i+1)+normps_diff(i+2)+normps_diff(i+3)+normps_diff(i+4)+normps_diff(i+5))/6;
        end
        normps_diff_samp=resample(normps_diff_smooth, length(freq_dqe), length(normps_diff_smooth), 0)';
        flatl=ceil(length(normps_diff_samp)*nps_squash);
        normps_diff_samp(1: flatl)=mean(normps_diff_samp(flatl+1:flatl+ceil(flatl/3)));
        nps_diff_slm = slmengine(freq_dqe, normps_diff_samp, 'plot', 'on', 'knots', 8);
        nps_diff=slmeval(freq_dqe, nps_diff_slm, 0)';
        figure;
        plot(normps_diff_x, log(normps_diff), 'b')
        hold on
        plot(freq_dqe, log(nps_diff), 'm')
        plot(freq_dqe, log(normps_red), 'r') 
        title('noise difference method')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Find DQE(0) and DQE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%PLOT the noise as a function of binning to find DQE(0)
%         if Find_DQE0 == 'y'
%         full_noise_area=map(1:1024, 1:1024);
%         %full_noise_area1=map(1:1024, 1:1024);
%         %full_noise_area=detrend(single(full_noise_area1));
%         s2_n2=zeros(10, 2);
%         var_between=zeros(10, 2);
%         var_within=zeros(10, 2);
%         
%         for i=1:5
%             stack_noise_squares=zeros(1024/(2^i), 1024/(2^i), 2^(i+1));
%             mean_vals=zeros(2^(i)^2, 1);
%             var_within_pre=zeros(2^(i)^2, 1);
%             indices=zeros(2^(i)^2, 2);
%             mean_square=zeros(2^i);
%             for j=1:2^(i)^2
%                  rindex=floor((j-1)/(2^i));
%                  cindex=mod((j-1), (2^i));
%                  sqsize=size(stack_noise_squares, 1);
%                 stack_noise_squares(:, :, j)=full_noise_area(rindex*sqsize+1:(rindex+1)*sqsize, ...
%                     cindex*sqsize+1:(cindex+1)*sqsize);
%                 tmp_square=stack_noise_squares(:,:,j);
%                 mean_vals(j)=mean(tmp_square(:));
%                 var_within_pre(j)=var(tmp_square(:));
%                 indices(j, 1)=rindex;
%                 indices(j, 2)=cindex;
%             end          
%             var_between(i, 1)=var(mean_vals(:))*(1024/(2^i))^2;
%             var_between(i, 2)=(1024/(2^i))^2;
%             if i==3
%                 for j=1:2^(i)^2
%                     mean_square(indices(j, 2)+1, indices(j, 1)+1)=mean_vals(j);
%                 end
%                 figure;
%                 image(mean_square, 'CDataMapping', 'scaled')
%                 colormap(gray)
%                 title('noise detrending?')
%             end
%             var_within(i, 1)=mean(var_within_pre(:));
%             var_within(i,2)=1024/(2^i);
%             s2_n2(i, 1)=(var_between(i));
%             %goes from 512 to 1
%             s2_n2(i, 2)=(1024/(2^i));
%         end
%         
%         s2_n2_copy=s2_n2;
%         
%         for i=1:3
%             counts=1;
%             sqside=1024/(2^i);
%             mean_vals=zeros(2^(i)^2, 1);
%             for j=sqside/2+1:4:1024-sqside/2
%                 for k=sqside/2+1:4:1024-sqside/2
%                     n_square=full_noise_area(j-sqside/2:j+sqside/2, k-sqside/2:k+sqside/2);
%                     mean_vals(counts)=mean(n_square(:));
%                     counts=counts+1;
%                 end
%             end
%             var_between(i, 1)=var(mean_vals(:))*(1024/(2^i))^2;
%             var_between(i, 2)=(1024/(2^i))^2;
%             s2_n2(i, 1)=var_between(i);
%         end         
%         
%         figure;
%      
%         plot(s2_n2(:, 2), s2_n2(:, 1), 'b')
%         hold on
%         plot(s2_n2_copy(:, 2), s2_n2_copy(:, 1), 'g')
%         hold off
%         title('S_n^2/n^2')
%         end
        %"For an incident beam obeying Poisson statistics, with on average n
        %electrons per pixel and giving an average output signal of d_n (so that
        %the average signal per primary electron is d_n/n)
        %DQE(w)=(d_n^2 * MTF(w)^2)/(n * NPS(w))"
        %McMullan et al 2009
        %Following derivation, DQE(w)=(DQE(0)*MTF^2)/(N~)
        %where DQE(0)=d_n^2/(n*NPS(0))
        %and N~ =NPS(w)/NPS(0)
        %Dimensional analysis: d_n/n is unitless; d_n is 1/pixel
%         if Find_DQE0 == 'y'
%             %minimap=double(map(1:1024, 1:1024));
%             maxshift=20;
%             %First find baseline variance
%             vari_base=0;
%             xc=0;
%             std1=0;
%             std2=0;
%             for rshift=0:maxshift-1
%                 for cshift=1:maxshift-1
%                     if rshift == maxshift-1 || cshift ==maxshift-1
%                     for i=1:size(map_meanvals, 1)-rshift
%                         for j=1:size(map_meanvals, 2)-cshift
%                             p1=map_meanvals(i, j);
%                             p2=map_meanvals(i+rshift, j+cshift);
%                             if p1 ~= noise_mean && p2 ~= noise_mean
%                             xc=xc+(p1-noise_mean)*(p2-noise_mean);
%                             std1=std1+(p1-noise_mean)^2;
%                             std2=std2+(p2-noise_mean)^2;
%                             end
%                         end                    
%                     end
%                     xc=xc/(sqrt((std1*std2)));
%                     vari_base=vari_base+xc; 
%                     end  
%                 end
%             end
%             vari_base=vari_base/(2*maxshift-1);
%             %Take minimap and shift a copy of it in increments from 1 to 20
%             %in two directions, calculating the cross corr at each shift.
%             xc=0;
%             vari=0;
%             std1=0;
%             std2=0;
%             for rshift=0:maxshift-1
%                 for cshift=1:maxshift-1
%                     for i=1:size(map_meanvals, 1)-rshift
%                         for j=1:size(map_meanvals, 2)-cshift
%                             if p1 ~= noise_mean && p2 ~= noise_mean
%                             p1=map_meanvals(i, j);
%                             p2=map_meanvals(i+rshift, j+cshift);
%                             xc=xc+(p1-noise_mean)*(p2-noise_mean);
%                             std1=std1+(p1-noise_mean)^2;
%                             std2=std2+(p2-noise_mean)^2;
%                             end
%                         end
%                     end
%                     xc=xc/(sqrt((std1*std2)));
%                     %exclude long-range correlations
%                     if xc >= 3*vari_base
%                     vari=vari+xc;
%                     end
%                 end
%             end
%             noise_rms=noise_rms*sqrt(1+8*vari)
%             vari_base
%             if vari_base >= 5/length(map_meanvals)
%                 disp('Warning--possible existence of long-range correlations in background from faulty gain correction or illumination gradient')
%             end
% 
%         %DQE0=noise_mean*gain_fac/mean(s2_n2(1:4))
%         DQE0=noise_mean*gain_fac/noise_rms^2;
%         end
 %sometimes the sizes are inexplicably off by one. This fixes that.
        if length(mtf_fh_samp) ~= length(normps_red)
            mtf_fh_temp=mtf_fh_samp;
            mtf_temp=mtf_samp;
            mtf_fh_samp=zeros(size(normps_red));
            mtf_samp=zeros(size(normps_red));
%             disp('in here')
%             size(mtf_fh_samp)
%             size(mtf_samp)
%             size(normps_red)
            for i=1:length(normps_red)
                mtf_fh_samp(i)=mtf_fh_temp(i);
                mtf_samp(i)=mtf_temp(i);
            end
        end
        DQE=(mtf_samp.^2)*noise_mean*gain_fac./normps_red;
        DQE_slm=(mtf.^2)*noise_mean*gain_fac./normps_red;
        if use_noise_file == 'y'
        size(nps_diff)
        size(mtf)
            DQE_diff=(mtf.^2)*noise_mean*gain_fac./nps_diff;
        end
        DQE_filt1=ifft(fftshift(fftshift(fft(DQE)).*hamming(length(DQE))));
        DQE_filt=sqrt(real(DQE_filt1).^2+imag(DQE_filt1).^2);
        DQE_filt(1)=DQE_filt(2);
        DQE_filt(length(DQE_filt))=DQE_filt(length(DQE_filt)-1);
%         DQE_ham=(mtf_fh_samp.^2)*noise_mean*gain_fac./normps_red;
%         DQE_ham_filt1=ifft(fftshift(fftshift(fft(DQE_ham)).*hamming(length(DQE_ham))));
%         DQE_ham_filt=sqrt(real(DQE_ham_filt1).^2+imag(DQE_ham_filt1).^2);
%         DQE_ham_filt(1)=DQE_ham_filt(2);
%         DQE_ham_filt(length(DQE_ham_filt))=DQE_ham_filt(length(DQE_ham_filt)-1);

        output_DQE=zeros(length(DQE_filt), 2);
        output_DQE(:, 1)=(0:1:length(DQE_filt)-1)/(length(DQE_filt)-1);
        output_DQE(:, 2)=DQE_filt;
        if use_noise_file == 'y'
            output_DQE(:, 2)=DQE_diff;
        end
        figure;
        plot(freq_dqe, DQE);
        hold on
        %plot((0:1:length(DQE_ham)-1)/(length(DQE_ham)-1), DQE_ham, 'g');
        plot(freq_dqe, DQE_filt, 'c');
        plot(freq_dqe, DQE_slm, 'r');
        if use_noise_file == 'y'
            plot(freq_dqe, DQE_diff, 'm');
        end
        %plot((0:1:length(DQE_ham_filt)-1)/(length(DQE_ham_filt)-1), DQE_ham_filt, 'r');
        axis([0 1 0 1])
        %scatter(0, DQE0)
        title('DQE blue/cyan is DQE from unsmoothed LSF. Red is from smoothed LSF and MTF. Magenta is smoothed and from noise diff. ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Find DQE MEthod 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 2
        %Mode 2: calculate MTF from FT of envelope.
        N=length(map);
        switch camera_type
            case 'f'
                map=map/1000;
                map=double(map);
            case {'u', 'o'}
                %map=map/10;
                map=double(map);
            otherwise
        end
        %Binarize and mask
        av_val=mean(mean(map(:)));
        std_val=std(single(map(:)));
        binmap=zeros(N);
        disp('found mean')
        background_acc=double([0, 0]);
        pointer_acc=double([0, 0]);
        for i=1:N^2
            if (map(i) >= (av_val-std_val))
               binmap(i)=1;
               background_acc(1)=background_acc(1)+(map(i));
               background_acc(2)=background_acc(2)+1;
            else
                binmap(i)=0;
                pointer_acc(1)=pointer_acc(1)+(map(i));
                pointer_acc(2)=pointer_acc(2)+1;
            end
        end
        background_av=background_acc(1)/background_acc(2);
        pointer_av=pointer_acc(1)/pointer_acc(2);
        disp('binarized')
        image(binmap', 'CDataMapping', 'scaled');
        caxis([0, 1])
        colormap(gray)
        axis equal
        set(gca, 'YDir', 'normal')
%         if MRC_write == 'y'
%         WriteMRC(binmap, 1, 'binary.mrc');
%         end
        figure;
        diffmat1=abs(diff(binmap, 1, 1));
        diffmat2=abs(diff(binmap, 1, 2));
        diffmat=zeros(N-1);
        for i=1:(N-1)
            for j=1:(N-1)
                if diffmat1(i, j)==1 || diffmat2(i, j)==1
                    diffmat(i, j)=1;
                end
            end
        end
        image(diffmat', 'CDataMapping', 'scaled');
        caxis([0, 1])
        set(gca, 'YDir', 'normal')
        axis equal
        colormap(gray)
%         if MRC_write == 'y'
%         WriteMRC(diffmat, 1, 'derivmap.mrc');
%         end
        %create mask by making a square around every point in the diffmat
        %if a point has less than five other points around it, it's no doubt an outlier
        %and should be thrown out.
        mask=zeros(N);
        mrad=20; %mask radius
        cew=40; %cosine edge width
        cos_x=(0:(pi/cew):pi);
        cos_edge=(1+cos(cos_x))/2;
        
        for i=1:N-1
            for j=1:N-1
                if diffmat(i, j) == 1
                    isone=0;
                    for k=-mrad:mrad
                    if (i+k <= N-1 && i+k >= 1)
                    for l=-mrad:mrad
                    if (j+l <= N-1 && j+l >= 1)
                        if diffmat(i+k, j+l) == 1 && (k ~=0 || l ~= 0)
                            isone=isone+1;
                        end
                    end; end; end; end
                    if isone <= 6
                        diffmat(i, j)=0;
                    else
                        for k=-mrad-cew:mrad+cew
                        %ifs to keep from checking out of bounds    
                        if (i+k <= N-1 && i+k >= 1)
                        for l=-mrad-cew:mrad+cew
                        if (j+l <= N-1 && j+l >= 1)
                        %mask is circle inscribed in the square
                        if (k^2+l^2<=mrad^2)
                            mask(i+k, j+l) = 1;
                        %cos edge extends outwards. Don't overwrite high
                        %mask vals with lower ones
                        elseif (k^2+l^2 > mrad^2 && k^2+l^2<=(mrad+cew)^2)
                            if mask(i+k, j+l) ~= 1 && ...
                                    cos_edge(ceil(sqrt(k^2+l^2))-mrad) > mask(i+k, j+l)
                            mask(i+k, j+l)=cos_edge(ceil(sqrt(k^2+l^2))-mrad);
                            end    
                        end;end;end;end;end
                    end
                end
            end
        end
        figure;
        image(mask', 'CDataMapping', 'scaled');
        caxis([0, 1])
        set(gca, 'YDir', 'normal')
        colormap(gray)
        axis equal
%         if MRC_write == 'y'
%         WriteMRC(mask, 1, 'mask.mrc');
%         end
        %Apply mask
        masked_im=zeros(N);
        for i=1:N-1
            for j=1:N-1
                if (mask(i, j) > 1 || mask(i, j) < 0)
                    disp('unexpected value in mask')
                    mask(i, j);
                end
                %if mask(i, j)==1
                %    masked_im(i, j)=map(i, j);
                %else
                    if binmap(i, j)==1
                        masked_im(i, j)=mask(i, j)*map(i, j)+(1-mask(i, j))*background_av;
                    else
                        masked_im(i, j)=mask(i, j)*map(i, j)+(1-mask(i, j))*pointer_av;
                    end
                %end
            end
        end
        for i=1:N
            if mask(i, N) == 1
                masked_im(i, N)=background_av;
            else
                masked_im(i, N)=pointer_av;
            end
            if mask(N, i) == 1
                masked_im(N, i)=background_av;
            else
                masked_im(N, i)=pointer_av;
            end
        end
        %Now for comparison purposes make the mask have the average values
        %instead of 1 and 0
        binpw=zeros(N);
        for i=1:N*N
            if binmap(i)==1
                binpw(i)=background_av;
            else
                binpw(i)=pointer_av;
            end
        end
        figure;
        image(masked_im', 'CDataMapping', 'scaled');
        colormap(gray)
        set(gca, 'YDir', 'normal')
        axis equal
        pwspec_mask=(abs(fftshift(fft2(masked_im))).^2)/(N^2);
        
        figure;
        subplot(1, 3, 1), image(log(pwspec_mask)', 'CDataMapping', 'scaled');
        colormap(gray)
        title('from map with perimeter-shaped mask')
        set(gca, 'YDir', 'normal')
        axis equal
%         if MRC_write == 'y'
%         WriteMRC(log(pwspec_mask), 1, 'log_pw_mask.mrc')
%         end
        pwspec_bin=(abs(fftshift(fft2(binpw))).^2)/(N^2);
        subplot(1, 3, 2), image(log(pwspec_bin)', 'CdataMapping', 'scaled');
        title('from binary map')
        set(gca, 'YDir', 'normal')
        axis equal
%         if MRC_write == 'y'
%         WriteMRC(log(pwspec_bin), 1, 'log_pw_bin.mrc')
%         end
        pwspec_raw=(abs(fftshift(fft2(map))).^2)/(N^2);
        subplot(1, 3, 3), image(log(pwspec_raw)', 'CdataMapping', 'scaled');
        title('from raw map')
        set(gca, 'YDir', 'normal')
        axis equal
%         if MRC_write == 'y'
%         WriteMRC(log(pwspec_raw), 1, 'log_pw_raw.mrc')
%         end
        
%         %create cross-shaped mask, 100 pixels wide
%         %mask_width=100;
%         mask_width=ceil(N/40);
%         %cross_mask=zeros(N);
%         cross_mask=ones(N);
%         cross_mask(N/2+1-mask_width/2:N/2+1+mask_width/2, :)=0;
%         cross_mask(:, N/2+1-mask_width/2:N/2+1+mask_width/2)=0;
%         %cross_mask(N/2+1-mask_width/2:N/2+1+mask_width/2, :)=1;
%         cross_mask(N/2+1, :)=0;
        
        %Back to the actual data
        %apply cross-shaped mask
%         pwspec_mask2=pwspec_mask.*cross_mask;
%         pwspec_bin2=pwspec_bin.*cross_mask;
% 
%         %Find mean of counts not behind mask:
%         mean_acc=0;
%         counts=0;
% %  q       exp_mask=zeros(N);
%         for i=1:(N*N)
%             if pwspec_bin2(i) ~= 0
%                 mean_acc=mean_acc+pwspec_bin2(i);
%                 counts=counts+1;
%             end
%         end
%         
%         av_val=mean_acc/counts
        %hist_x=min(pwspec_mask2(:)):(max(pwspec_mask2(:))-min(pwspec_mask2(:)))/100:max(pwspec_mask2(:));
        %hist_x
%         [n, xout]=hist(pwspec_mask2(:), 100); 
%         figure;
%         imhist(pwspec_bin2)
%         
%         for i=1:(N*N)
%             if pwspec_bin2(i) >= av_val*0.05
%                 exp_mask(i)=1;
%             end
%         end
        
%         figure;
%         image(exp_mask', 'CdataMapping', 'scaled');
%         title('exp_mask')
%         set(gca, 'YDir', 'normal')
%         axis equal
%         downs_mask=imresize(exp_mask, 1/200);
%         for i=1:(size(downs_mask, 1)*size(downs_mask, 2))
%             if downs_mask(i) >= 0.5
%                 downs_mask(i)= 1;
%             else
%                 downs_mask(i)=0;
%             end
%         end
%         new_mask=imresize(downs_mask, [N, N]);
%         figure;
%         image(downs_mask', 'CDataMapping', 'scaled');
%         title('small fourier mask (downs_map)')
%         set(gca, 'YDir', 'normal')
%         colormap(gray)
%         axis equal
%         caxis([0 1]);
%         figure;
%         image(new_mask', 'CDataMapping', 'scaled');
%         title('new fourier mask (new_map)')
%         set(gca, 'YDir', 'normal')
%         colormap(gray)
%         axis equal
%         caxis([0 1]);
        
%         new_mask_bw=im2bw(new_mask, 0.5);
%         pwspec_mask3=pwspec_mask2.*new_mask_bw;
%         pwspec_bin3=pwspec_bin2.*new_mask_bw;
%         %TEMPORARY
%         pwspec_mask3=pwspec_mask2;
%         pwspec_bin3=pwspec_bin2;
        
%         figure;
%         subplot(1, 3, 1), image(log(pwspec_bin3)', 'CdataMapping', 'scaled');
%         title('from binary map')
%         set(gca, 'YDir', 'normal')
%         axis equal
%         if MRC_write == 'y'
%         WriteMRC(log(pwspec_bin3), 1, 'cross_ps_bin.mrc');
%         end
%         subplot(1, 3, 2), image(log(pwspec_mask3)', 'CDataMapping', 'scaled');
%         title('from map with perimeter-shaped mask')
%         set(gca, 'YDir', 'normal')
%         axis equal
%         if MRC_write == 'y'
%         WriteMRC(log(pwspec_mask3), 1, 'cross_ps_mask.mrc');
%         end
%         subplot(1, 3, 3), image(cross_mask', 'CDataMapping', 'scaled')
%         set(gca, 'YDir', 'normal')
%         axis equal
%         caxis([0 1])
%         colormap(gray)
        
        %Radial average and plot
        counts_mask=zeros((N/2)+1, 1);
        ps_mask=zeros((N/2)+1, 1);
        ps_mask_weights=zeros((N/2)+1, 1);
        counts_bin=zeros((N/2)+1, 1);
        ps_bin=zeros((N/2)+1, 1);
        ps_bin_weights=zeros((N/2)+1, 1);
        counts_raw=zeros((N/2)+1, 1);
        ps_raw=zeros((N/2)+1, 1);
        
        for i=1:N
            for j=1:N/2+1
                radius=ceil(sqrt((i-(N/2+1))^2+(j-(N/2+1))^2));
                if radius <= length(ps_mask)-1
                    if pwspec_bin(i, j) ~= 0
                    counts_mask(radius+1)=counts_mask(radius+1)+1;
                    ps_mask(radius+1)=ps_mask(radius+1)+pwspec_mask(i, j);
                    ps_mask_weights(radius+1)=ps_mask_weights(radius+1) ...
                        +pwspec_mask(i, j)*pwspec_bin(i,j)^5;
                    counts_bin(radius+1)=counts_bin(radius+1)+1;
                    ps_bin(radius+1)=ps_bin(radius+1)+pwspec_bin(i, j);
                    ps_bin_weights(radius+1)=ps_bin_weights(radius+1) ...
                        +pwspec_bin(i, j)^6;
                    counts_raw(radius+1)=counts_raw(radius+1)+1;
                    ps_raw(radius+1)=ps_raw(radius+1)+pwspec_raw(i, j);
                    end
                end
            end
        end

        normps_mask=bsxfun(@rdivide, ps_mask, counts_mask);
        normps_bin=bsxfun(@rdivide, ps_bin, counts_bin);
        normps_raw=bsxfun(@rdivide, ps_raw, counts_raw);
        n_bin_weights=bsxfun(@rdivide, ps_bin_weights, counts_bin);
        n_mask_weights=bsxfun(@rdivide, ps_mask_weights, counts_mask);
        
        freq=(0:1:length(normps_mask)-1)/(length(normps_mask)-1);
        figure;
        ratio_ps=zeros(N/2+1, 1);
        ratio_ps_weights=zeros(N/2+1, 1);
        for i=1:N/2+1
            if normps_bin(i) ~= 0
                ratio_ps(i)=normps_mask(i)/normps_bin(i);
            end
        end
        for i=1:N/2+1
            if n_bin_weights(i) ~=0
                ratio_ps_weights(i)=n_mask_weights(i)/n_bin_weights(i);
            end
        end
        
        subplot(1, 2, 1), plot(freq, normps_mask, 'b', ...
            freq, normps_bin, 'r', freq, ratio_ps.^(0.5), 'g')
        title('masked=blue, binary=red, mtf=green')
        axis([-Inf Inf 0 1])
        subplot(1, 2, 2), plot(freq, n_mask_weights, 'b', ...
            freq, n_bin_weights, 'r', freq, ratio_ps_weights.^(0.5), 'g')
        title('MTF. masked=blue, binary=red, mtf=green. With weights.')
        axis([-Inf Inf 0 1])
        
        mtf=ratio_ps_weights.^(0.5);
        disp('size rps')
        size(ratio_ps_weights)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Noise Power Spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

%Divide image into npatches patches of size spatches, marking the center of
%each patch with an areadot
        spatches=512; %Must be even
        npatches=floor(size(map, 1)/spatches)*floor(size(map, 2)/spatches);
        noise_areas=zeros(spatches, spatches, npatches);
        area_dots=zeros(npatches, 2);
        noise_fts=zeros(spatches, spatches, npatches);
        padsq_length=length(fftshift(fft2(single(noise_areas(:,:, 1)), spatches*4, spatches*4)));
        noise_fts_pad=zeros(padsq_length, padsq_length, npatches);
        nps_stack=zeros(spatches, spatches, npatches);
        nps_stack_pad=zeros(padsq_length, padsq_length, npatches);
        noise_stats=zeros(npatches, 2);
        noise_include=ones(npatches, 1);
        index_convert=zeros(npatches, 2);
        
        %Create a cose-edge mask to go around all of the square patches.
        cew=10;
        cos_x=(0:(pi/cew):pi);
        cos_edge=1-((1+cos(cos_x))/2);
        cos_mask=ones(spatches);
        for i=1:cew
            for j=1:spatches
            if cos_mask(i, j) > cos_edge(i)
                cos_mask(i, j)=cos_edge(i);
                cos_mask(j, i)=cos_edge(i);
                cos_mask(spatches-i, j)=cos_edge(i);
                cos_mask(j, spatches-i)=cos_edge(i);
            end
        end
        end
            
           
        mask_weight=sum(sum(cos_mask.^2));
        size(map)
        for i=1:floor(size(map, 1)/spatches)
            for j=1:floor(size(map, 2)/spatches)
                index1d=(i-1)*floor(size(map, 1)/spatches)+j;
                index_convert(index1d, 1)=j;
                index_convert(index1d, 2)=i;
                noise_areas(:,:,index1d)=map((i-1)*spatches+1:(i-1)*spatches+spatches, ...
                    (j-1)*spatches+1:(j-1)*spatches+spatches, 1);
                noise_mean_tmp=mean(mean(noise_areas(:, :, index1d)));
                noise_areas(:, :, index1d)=noise_areas(:, :, index1d).*cos_mask+noise_mean_tmp*(1-cos_mask);
                area_dots(index1d, 1)=spatches*(j-1)+spatches/2;                
                area_dots(index1d, 2)=spatches*(i-1)+spatches/2;
                noise_fts(:,:,index1d)=fftshift(fft2(single(noise_areas(:,:,index1d))));
              %  noise_fts_pad(:, :, index1d)=fftshift(fft2(single(noise_areas(:,:, index1d)), spatches*4, spatches*4));         
                nps_stack(:, :, index1d)=(abs(noise_fts(:,:,index1d)).^2)/(spatches*spatches);
               % nps_stack_pad(:, :, index1d)=(abs(noise_fts_pad(:,:,index1d)).^2)/(spatches^2);
                noise_tmp=noise_areas(:, :, index1d);
                noise_stats(index1d, 1)=mean(noise_tmp(:));
                noise_stats(index1d, 2)=var(noise_tmp(:));                
            end
        end

        %subplot(2, 2, 4), image(log(nps_stack(:, :, 3)), 'CDataMapping', 'scaled')
        %colormap(gray)

%Determine which patches contain pointer by finding patches with 
%significantly lower mean values. Do statistics on only largest half of values
%to avoid pointer. Mark pointer patches with a red dot and set noise_include to 0.
        noise_mean_sort=sort(noise_stats(:, 1), 'descend');
        noise_mean_mean=mean(noise_mean_sort(1:ceil(length(noise_mean_sort)/2)));
        noise_mean_stdev=std(noise_mean_sort(1:ceil(length(noise_mean_sort)/2)));
        red_dots=0;
        for i=1:npatches
            if noise_stats(i, 1)<=noise_mean_mean-3*noise_mean_stdev
                noise_include(i)=0;
                red_dots=red_dots+1;
            end
        end
        red_dot_coords=zeros(red_dots, 2);
        j=1;
        for i=1:npatches
            if noise_include(i) == 0
                red_dot_coords(j, 1) = area_dots(i, 1);
                red_dot_coords(j, 2) = area_dots(i, 2);
                j= j+1;
            end
        end    
        figure;
        image(mapt, 'CDataMapping', 'scaled')
        colormap(gray)
        axis equal
        set(gca, 'YDir', 'normal')
        title('red dots on the pointer, blue dots on noise patches')
        hold on
        scatter(area_dots(:, 1), area_dots(:, 2))
        scatter(red_dot_coords(:, 2), red_dot_coords(:, 1), 'r')
%Average NPS; average noise value
        tmpsum_ps=zeros(size(nps_stack(:,:,1)));
        tmpsum_real=zeros(size(noise_areas(:,:,1)));
        tmpsum_squared=zeros(size(noise_areas(:,:,1)));
        tmpcnts=0;
        for i=1:length(noise_include)
            if noise_include(i) == 1
                tmpsum_ps=tmpsum_ps+nps_stack(:,:,i);
%                   figure;
%                   image(log(nps_stack(:, :, i)), 'CDataMapping', 'scaled')
%                   colormap(gray)
                tmpsum_real=tmpsum_real+noise_areas(:,:,i);
                tmpsum_squared=tmpsum_squared+noise_areas(:,:,i).^2;
                tmpcnts=tmpcnts+1;
            end
        end
        nps=tmpsum_ps/tmpcnts;
        noise_squared=sum(tmpsum_squared(:))/tmpcnts/length(noise_areas)^2;
        noise_mean=mean(mean(tmpsum_real/tmpcnts))
        noise_rms=sqrt(noise_squared-noise_mean^2)
       
%Radial average and plot
% 
%         tmp=fftshift(fft2(single(map(2500:3524, 2500:3524))));
%         nps=(abs(tmp).^2)/(1024*1024);
%         figure;
%         image(log(nps), 'CDataMapping', 'scaled');
%         colormap(gray);
        
        enditer=spatches^2;
        cr=spatches/2+1;
        cc=spatches/2+1;
        counts=zeros(spatches/2+1, 1);
        ps=zeros(spatches/2+1, 1);
        for i=1:enditer
            [r c]=ind2sub([spatches, spatches], i);
            radius=ceil(sqrt((r-cr)^2+(c-cc)^2));
            if radius <= length(ps)-1 && radius > 0
                counts(radius)=counts(radius)+1;
                ps(radius)=ps(radius)+nps(r, c);
            end
        end


        normps=bsxfun(@rdivide, ps, counts);
        normps(length(normps))=normps(length(normps)-2);
        normps(length(normps)-1)=normps(length(normps)-2);
        
        normps_hammed=normps;
         ham_x=(0:1:length(normps_hammed)-1)/(length(normps_hammed)-1);

        flatl=ceil(length(normps_hammed)*nps_squash)
        normps_hammed(1: flatl)=mean(normps_hammed(flatl+1:flatl+ceil(flatl/3)));
      
        %Smooth the power spectrum
        normps_smooth=normps_hammed;
%         for i=1:length(normps_hammed)-6
%             normps_smooth(i)=(normps_hammed(i)+normps_smooth(i+1)+normps_smooth(i+2)+normps_smooth(i+3)+normps_smooth(1+4)+normps_smooth(i+5))/6;
%         end
        
        nps_slm = slmengine(ham_x, normps_smooth, 'plot', 'on');
        nps_slm_eval=slmeval(ham_x, nps_slm, 0);
        normps_filt=medfilt1(nps_slm_eval, 6);
        %normps_filt=nps_slm_eval;
        %resize nps by resampling so it is the same as the mtf
        if(length(normps_filt) >= length(mtf)) 
            disp('the resampling is going from')
            disp('checkpoint')
            length(normps_filt)
            length(mtf)
            normps_red=resample(normps_filt, length(xfreq), length(normps_filt), 0)';
        else
            disp('the resampling is going from')
            disp('checkpoint 2')
            size(mtf)
            size(normps_filt)
            x_samp=(0:1:length(normps_filt)-1)/(length(normps_filt)-1);
            x_resamp=(0:1:length(mtf)-1)/(length(mtf)-1);
            normps_red=interp1(x_samp, normps_filt, x_resamp)';
        end
        
        normps_red(1)=normps_red(2);
        freq_dqe=(0:1:length(normps_red)-1)/(length(normps_red)-1);
        freq_dqe_long=(0:1:length(normps_filt)-1)/length(normps_filt);
        freq_dqe_smooth=(0:1:length(normps_smooth)-1)/length(normps_smooth);
        
        figure;
        plot(freq_dqe_long, log(normps));
        hold on
        title('log(Noise Power Spectrum)')
        plot(freq_dqe_long, log(normps_filt), 'r');
        %plot(ham_x, log(normps_hammed), 'g');
        plot(freq_dqe, log(normps_red), 'or')
        plot(freq_dqe_smooth, log(normps_smooth), 'og')
        figure;
        plot(freq_dqe_smooth, normps_smooth/(normps_smooth(3)))
        title('NPS normalized to 1')
        
        %sometimes the sizes are inexplicably off by one. This fixes that.
        if length(mtf) ~= length(normps_red)
            mtf_temp=mtf;
            mtf=zeros(size(normps_red));
            disp('in here')
            for i=1:length(normps_red)
                mtf(i)=mtf_temp(i);
            end
        end

        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Find DQE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First, create an SSNR filter
        nps_x=(0:1:length(normps_bin)-1)'/(length(normps_bin)-1);
        flatl=0.01
        flatl2=ceil(flatl*(length(normps_bin)))
        normps_bin(1:flatl2)=mean(normps_bin(flatl2+1:flatl2*1.3));
        bin_slm = slmengine(nps_x, normps_bin, 'plot', 'on', 'knots', [0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1]);
        bin_slm_eval=slmeval(nps_x, bin_slm, 0);
        normps_mask(1:flatl2)=mean(normps_mask(flatl2+1:1.3*flatl2));
        mask_slm = slmengine(nps_x, normps_mask, 'plot', 'on', 'knots', [0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1]);
        mask_slm_eval=slmeval(nps_x, mask_slm, 0);
        normps_raw(1:flatl2)=mean(normps_raw(flatl2+1:1.3*flatl2));
        raw_slm = slmengine(nps_x, normps_raw, 'plot', 'on', 'knots', [0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1]);
        raw_slm_eval=slmeval(nps_x, raw_slm, 0);
          
        figure;
        plot(bin_slm_eval)
        hold on
        plot(mask_slm_eval, 'r');
        hold on
        plot(raw_slm_eval, 'g');

        normps_noise=raw_slm_eval-bin_slm_eval;
        for i=1:length(normps_noise)
            if normps_noise(i) < 0
                normps_noise(i)=0;
            end
        end
        SSNR=normps_bin./normps_noise;
        SSNR_filt=SSNR./(1+SSNR);
        for i=1:length(SSNR_filt)
            if isnan(SSNR_filt(i))
                SSNR_filt(i)=1;
            end
        end
    SSNR_filt
        size(SSNR_filt)
        size(normps_red)
        DQE=(mtf.^2)*noise_mean*gain_fac./normps_red;
        DQE_filt1=ifft(fftshift(fftshift(fft(DQE)).*hamming(length(DQE))));
        DQE_filt=sqrt(real(DQE_filt1).^2+imag(DQE_filt1).^2);
        DQE_filt(1)=DQE_filt(2);
        DQE_filt(length(DQE_filt))=DQE_filt(length(DQE_filt)-1);
        output_DQE=zeros(length(DQE_filt), 2);
        output_DQE(:, 1)=(0:1:length(DQE_filt)-1)/(length(DQE_filt)-1);
        output_DQE(:, 2)=DQE_filt;
        figure;
        plot(freq_dqe, DQE);
        hold on
        plot((0:1:length(DQE_filt)-1)/(length(DQE_filt)-1), DQE_filt, 'g');
        plot(freq_dqe, DQE.*SSNR_filt, 'r')
        axis([0 1 0 1])
        title('DQE')

end

end







