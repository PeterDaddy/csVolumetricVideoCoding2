clear;
close all;
clc;
% profile on
UpFolder = fileparts(pwd);
addpath(fullfile(UpFolder, 'libraries\l1magic\Optimization'));
addpath(fullfile(UpFolder, 'libraries\l1magich\Measurements'));
addpath(fullfile(UpFolder, 'libraries\l1magic\Data'));
addpath(fullfile(UpFolder, 'libraries\spgl1-1.9'));
addpath(fullfile(UpFolder, 'libraries\SL0'));
addpath(fullfile(UpFolder, 'libraries\TwIST_v2'));
addpath(fullfile(UpFolder, 'libraries\NESTA_v1.1'));
addpath(fullfile(UpFolder, 'libraries\L1_homotopy_v2.0'));
addpath(fullfile(UpFolder, 'libraries\lzw'));
addpath(fullfile(UpFolder, 'libraries\pureAC'));
addpath(fullfile(UpFolder, 'libraries\'));
addpath(fullfile(UpFolder, 'libraries/BallLabsAlgo'));
addpath(fullfile(UpFolder, 'sequences\'));
imageOriginalPath = 'C:\Users\jiray\OneDrive\Desktop\Doctor\Research\sequences\';
imageFiles = [dir(fullfile(imageOriginalPath,'*png'));
              dir(fullfile(imageOriginalPath,'*tiff'));
              dir(fullfile(imageOriginalPath,'*tif'));
              dir(fullfile(imageOriginalPath,'*gif'));
              dir(fullfile(imageOriginalPath,'*jpg'));
              dir(fullfile(imageOriginalPath,'*bmp'));
              dir(fullfile(imageOriginalPath,'*mat'))];
numFiles = length(imageFiles);
%___SIMULATION SETUPS___
simulation_parameter.macro_block_size                = 16;
simulation_parameter.sliding_window_size             = 4;
simulation_parameter.search_area                     = 8;
simulation_parameter.n                               = simulation_parameter.macro_block_size^2; % NOTE: small error still present after increasing m;
simulation_parameter.measurement_matrix_lists        = [simulation_parameter.macro_block_size^2*1];
simulation_parameter.measurement_matrix_construction = 'binary_walsh_zigzag';
simulation_parameter.reconstruction_algorithm        = 'l1qc_logbarrier';
simulation_parameter.transformation_algorithm        = 'ifwht';
simulation_parameter.color_mode                      = 'gray';

intra_inter_prediction_obj = intra_inter_prediction;
for matrix_depth = [64]
    matrix_depth
    simulation_parameter.m = matrix_depth;
    switch simulation_parameter.measurement_matrix_construction
        case 'binary_toeplitz'
            temp_matrix              = load('toeplitz_matrix.mat');
            simulation_parameter.phi = temp_matrix.ans(1:simulation_parameter.m,1:simulation_parameter.n);
        case 'binary_parlay'
            q = log2(simulation_parameter.n);
            if  sum(ismember(char(cellstr(num2str(q))),'.'))~=0
                disp('           Warning!...               ');
                disp('The size of Vector  must be in the shape of 2^N ..');
                return
            else
                for u = 1:simulation_parameter.n
                    binu = dec2bin(u-1,q);
                    for v = 1:simulation_parameter.n
                        binv = dec2bin(v-1,q);
                        temp = 0;
                        for i = 1:q
                            temp= temp + bin2dec(binu(i))*bin2dec(binv(q+1-i));
                        end
                       D(u,v)=(-1)^temp;
                    end
                end
            end
            simulation_parameter.phi = max(D(1:simulation_parameter.m,1:simulation_parameter.n),0);
        case 'binary_random'
            temp_matrix              = load('random_matrix.mat');
            simulation_parameter.phi = temp_matrix.simulation_parameter.phi(1:simulation_parameter.m,1:simulation_parameter.n);
        case 'binary_hadamard'
            hadamard_matrix          = hadamard(simulation_parameter.n);  
            simulation_parameter.phi = max(hadamard_matrix,0);
            sub_phi.x                = zeros(1, size(simulation_parameter.phi,1)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
            sub_phi.y                = zeros(1, size(simulation_parameter.phi,2)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
            simulation_parameter.phi = max(simulation_parameter.phi(1:simulation_parameter.m,1:simulation_parameter.n),0);
        case 'binary_walsh'
            hadamard_matrix          = hadamard(simulation_parameter.n);
            HadIdx                   = 0:simulation_parameter.n-1;       % Hadamard index
            M                        = log2(simulation_parameter.n)+1;    % Number of bits to represent the index
            binHadIdx                = fliplr(dec2bin(HadIdx,M))-'0';     % Bit reversing of the binary index
            binSeqIdx                = zeros(simulation_parameter.n,M-1); % Pre-allocate memory
            for k = M:-1:2
                % Binary sequency index
                binSeqIdx(:,k)       = xor(binHadIdx(:,k),binHadIdx(:,k-1));
            end
            SeqIdx                   = binSeqIdx*pow2((M-1:-1:0)');   % Binary to integer sequency index
            walshMatrix              = hadamard_matrix(SeqIdx+1,:);   % 1-based indexing
            simulation_parameter.phi = max(walshMatrix(1:simulation_parameter.m,1:simulation_parameter.n), 0);
        case 'binary_walsh_zigzag'
            hadamard_matrix          = hadamard(simulation_parameter.n);
            HadIdx                   = 0:simulation_parameter.n-1;        % Hadamard index
            M                        = log2(simulation_parameter.n)+1;    % Number of bits to represent the index
            binHadIdx                = fliplr(dec2bin(HadIdx,M))-'0';     % Bit reversing of the binary index
            binSeqIdx                = zeros(simulation_parameter.n,M-1); % Pre-allocate memory
            for k = M:-1:2
                % Binary sequency index
                binSeqIdx(:,k)       = xor(binHadIdx(:,k),binHadIdx(:,k-1));
            end
            SeqIdx                   = binSeqIdx*pow2((M-1:-1:0)');   % Binary to integer sequency index
            walshMatrix              = hadamard_matrix(SeqIdx+1,:);   % 1-based indexing
            simulation_parameter.phi = max(walshMatrix,0);
            sub_phi.x                = zeros(1, size(simulation_parameter.phi,1)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
            sub_phi.y                = zeros(1, size(simulation_parameter.phi,2)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
            simulation_parameter.phi = mat2cell(simulation_parameter.phi, sub_phi.x, sub_phi.y);
            %___Re-order to AC DC
            t=0;
            l=size(simulation_parameter.phi);
            sum_=l(2)*l(1);  %calculating the M*N
            for d=2:sum_
                c=rem(d,2);  %checking whether even or odd
                for i=1:l(1)
                    for j=1:l(2)
                        if((i+j)==d)
                            t=t+1;
                            if(c==0)
                                simulation_parameter.phi_temp(t,:) = simulation_parameter.phi{j,d-j}(:);
                            else          
                                simulation_parameter.phi_temp(t,:) = simulation_parameter.phi{d-j,j}(:);
                            end
                         end    
                     end
                 end
            end
            simulation_parameter.phi = max(simulation_parameter.phi_temp(1:simulation_parameter.m,1:simulation_parameter.n),0);
    end

    mat = simulation_parameter.phi;                   % Your sample matrix
    [r, c] = size(mat);                               % Get the matrix size
    imagesc((1:c)+0.5, (1:r)+0.5, mat);               % Plot the image
    colormap(gray);                                   % Use a gray colormap
    axis equal                                        % Make axes grid sizes equal
    set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
             'XLim', [1 c+1], 'YLim', [1 r+1], ...
             'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
    %___THETA___
    %___NOTE: Avoid calculating Psi (nxn) directly to avoid memory issues___
     simulation_parameter.theta = zeros(simulation_parameter.m,simulation_parameter.n);
    for theta_loop = 1:simulation_parameter.n
        ek = zeros(1,simulation_parameter.n);
        ek(theta_loop) = 1;
        switch simulation_parameter.transformation_algorithm
            case 'ifwht'
                simulation_parameter.psi = ifwht(ek)';
        end
        simulation_parameter.theta(:,theta_loop) = simulation_parameter.phi*simulation_parameter.psi;
    end
    
    intra_buffer_encoder = cell(136, 240, 3);
    inter_buffer_encoder = cell(136, 240, 3);
    for kk = 1:3
        for ii = 1:136
            for jj = 1:240
                intra_buffer_encoder{ii,jj,kk} = zeros([simulation_parameter.sliding_window_size simulation_parameter.sliding_window_size]);
                inter_buffer_encoder{ii,jj,kk} = zeros([simulation_parameter.sliding_window_size simulation_parameter.sliding_window_size]);
            end
        end
    end
    
    for frame_number = 1:1
        frame_number
        %___LOAD IMAGE___
        load_frame = imread(imageFiles(frame_number).name);
        if(strcmp(simulation_parameter.color_mode,'rgb') || strcmp(simulation_parameter.color_mode,'RGB'))
            frame = load_frame(:,:,:);
%             frame = padarray(frame, [4 2]);
            plane = 3;
        elseif(strcmp(simulation_parameter.color_mode,'gray') || strcmp(simulation_parameter.color_mode,'GRAY'))
            frame = double(rgb2gray(load_frame));
%             frame = padarray(frame, [4 2]);
            %frame = frame(1:720, 1:1280, :);
            plane = 1;
        else
            frame = load_frame;
            plane = 1;
        end

        %___RESET STATE___
        sub_block.x = zeros(1, size(frame,1)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
        sub_block.y = zeros(1, size(frame,2)/simulation_parameter.macro_block_size) + simulation_parameter.macro_block_size;
        for k = 1:plane
            frame_temp(:,:,k) = mat2cell(frame(:,:,k), sub_block.x, sub_block.y);
        end
        
        %___THE RANDOM PROJECTION___
        disp('Random Projection...');
        for k = 1:plane
            for i = 1:size(frame,1)/simulation_parameter.macro_block_size
                for j = 1:size(frame,2)/simulation_parameter.macro_block_size
                   one_block_image(:,:,k) = reshape(frame_temp{i,j,k},simulation_parameter.macro_block_size^2,1);
                   y.measurement{i,j,k}   = (BCS_encoder(double(one_block_image(:,:,k)), simulation_parameter.phi)); %___Sampling
                end
            end
        end
        disp('Random Projection Done');
        disp('Measurement plane extraction...');
        y.norm_sub_frame = 0;
        for k = 1:plane
            for i = 1:size(frame,1)/simulation_parameter.macro_block_size
                for j = 1:size(frame,2)/simulation_parameter.macro_block_size
                    for z = 1:simulation_parameter.m
                        y.sub_frame(i,j,k,z) = y.measurement{i,j,k}(z);
                    end
                end
            end
        end
        
        for k = 1:plane
            y.norm_sub_frame = y.sub_frame(:,:,k,z);
        end
        imwrite(ind2rgb(im2uint8(mat2gray(y.norm_sub_frame)), gray(256)), strcat('layer', num2str(z), '.bmp'))
        
        disp('Done');
        
        disp('Prediction template generating...');
        for k = 1:plane
            %___Get one sub_frame___
            sub_frame                                   = y.norm_sub_frame;
            %___Special_case padding for 4K___
            sub_frame(136,:)                            = 0;
            sub_block.x                                 = zeros(1, size(sub_frame,1)/simulation_parameter.sliding_window_size) + simulation_parameter.sliding_window_size;
            sub_block.y                                 = zeros(1, size(sub_frame,2)/simulation_parameter.sliding_window_size) + simulation_parameter.sliding_window_size;
            sub_frame_temp                              = mat2cell(sub_frame, sub_block.x, sub_block.y);
            for i = 1:size(sub_frame_temp,1)
                for j = 1:size(sub_frame_temp,2)
                    inter_buffer_encoder                = intra_buffer_encoder;
%                     if(frame_number == 1)
                        prediction_candidate            = intra_inter_prediction_obj.intra_prediction(sub_frame_temp, intra_buffer_encoder, simulation_parameter.sliding_window_size, i, j);
                        internal_residual{i,j,k}        = sub_frame_temp{i,j} - prediction_candidate;
                        prediction_template{i,j,k}      = prediction_candidate;
                        qp                                      = 16;
                        q_mtx                                   = [2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)); 
                                                                   2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)); 
                                                                   2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp));
                                                                   2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp))];
                        %___Do more about quantization___
                        quantization_coded_sub_frame{i,j,k}   = round(internal_residual{i,j,k}./q_mtx);
                        dequantization_coded_sub_frame{i,j,k} = quantization_coded_sub_frame{i,j,k}.*q_mtx;
                        intra_buffer_encoder{i,j,k}           = dequantization_coded_sub_frame{i,j,k} + prediction_candidate;
%                     else
%                         intra_prediction_candidate              = intra_inter_prediction_obj.intra_prediction(sub_frame_temp, intra_buffer_encoder, simulation_parameter.sliding_window_size, i, j);
%                         inter_prediction_candidate              = intra_inter_prediction_obj.inter_prediction(sub_frame_temp, inter_buffer_encoder, simulation_parameter.sliding_window_size, i, j);
%                         if(sum(sum(abs(sub_frame_temp{i,j}-intra_prediction_candidate)))<sum(sum(abs(sub_frame_temp{i,j}-inter_prediction_candidate))))
%                             prediction_candidate                = intra_prediction_candidate;
%                             residual_sub_frame{i,j,k,z}         = sub_frame_temp{i,j} - prediction_candidate;
%                         else
%                             prediction_candidate                = inter_prediction_candidate;
%                             residual_sub_frame{i,j,k,z}         = sub_frame_temp{i,j} - prediction_candidate;
%                         end
%                         %___Generate DCT 2D Matrix___
%                         dct2dmx                                 = dctmtx(4);
%                         dct_coded_sub_frame{i,j,k,z}            = dct2dmx*residual_sub_frame{i,j,k,z}*dct2dmx';
%                         q_mtx                                   = [1  2  4  8; 
%                                                                    2  4  8 16; 
%                                                                    4  8 16 32;
%                                                                    8 16 32 64];
%                         %___Do more about quantization___
%                         quantization_coded_sub_frame{i,j,k,z}   = round(dct_coded_sub_frame{i,j,k,z}./q_mtx);
%                         dequantization_coded_sub_frame{i,j,k,z} = quantization_coded_sub_frame{i,j,k,z}.*q_mtx;
%                         invdct_coded_sub_frame{i,j,k,z}         = dct2dmx'*dequantization_coded_sub_frame{i,j,k,z}*dct2dmx;
%                         coded_sub_frame{i,j,k,z}                = invdct_coded_sub_frame{i,j,k,z} + prediction_candidate;
%                         %__Next iteration__
%                         %___Big Bug___it should not be cell but matrix
%                         %in set of z
%                         intra_buffer_encoder{i,j,k}             = invdct_coded_sub_frame{i,j,k,z} + prediction_candidate;
%                     end
                end
            end
        end
        disp('Done');

        disp('Compressing...');
        for k = 1:plane
            for z = 1:simulation_parameter.m
                %___Get one sub_frame___
                sub_frame                   = y.sub_frame(:,:,k,z);
                %___Special_case padding for 4K___
                sub_frame(136,:)            = 0;
                sub_block.x                 = zeros(1, size(sub_frame,1)/simulation_parameter.sliding_window_size) + simulation_parameter.sliding_window_size;
                sub_block.y                 = zeros(1, size(sub_frame,2)/simulation_parameter.sliding_window_size) + simulation_parameter.sliding_window_size;
                sub_frame_temp              = mat2cell(sub_frame, sub_block.x, sub_block.y);
                for i = 1:size(sub_frame_temp,1)
                    for j = 1:size(sub_frame_temp,2)
                        residual_sub_frame{i,j,k,z}             = sub_frame_temp{i,j} - prediction_template{i,j,k};
                        qp                                      = 16;
                        q_mtx                                   = [2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)); 
                                                                   2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)); 
                                                                   2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp));
                                                                   2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp)) 2^(log2(qp))];
                        %___Do more about quantization___
                        quantization_coded_sub_frame{i,j,k,z}   = round(residual_sub_frame{i,j,k,z}./q_mtx);
                        dequantization_coded_sub_frame{i,j,k,z} = quantization_coded_sub_frame{i,j,k,z}.*q_mtx;
                        coded_sub_frame{i,j,k,z}                = dequantization_coded_sub_frame{i,j,k,z} + prediction_template{i,j,k};
                    end
                end
            end 
        end
        disp('Done');
        
        disp('Entropy measuring...');
        translation_temp = cell2mat(quantization_coded_sub_frame);
        %___Trancate to original size
        trancate_translation_temp = translation_temp(1:size(frame,1)/simulation_parameter.macro_block_size, ...
                                                     1:size(frame,2)/simulation_parameter.macro_block_size, ...
                                                     1:plane, ...
                                                     1:simulation_parameter.m);
        entropy_bpp_compressed = 0;
        for k = 1:plane
            for i = 1:size(trancate_translation_temp,1)
                for j = 1:size(trancate_translation_temp,2)
                    for z = 1:simulation_parameter.m
                        y.element_plane_translation{i,j,k}(z,1) = trancate_translation_temp(i,j,k,z);
                    end
                    if(max(y.element_plane_translation{i,j,k}) == min(y.element_plane_translation{i,j,k}))
                        entropy_bpp_compressed = entropy_bpp_compressed + (1/(size(frame,1)*size(frame,2)));
                    else
                        entropy_bpp_compressed = entropy_bpp_compressed + (Measurement_Entropy(y.element_plane_translation{i,j,k},(size(frame,1)*size(frame,2))));
                    end
                end
            end
        end
        disp('Done');
        
        disp('Element plane to vector translating...');
        translation_temp = cell2mat(coded_sub_frame);
        %___Trancate to original size
        trancate_translation_temp = translation_temp(1:size(frame,1)/simulation_parameter.macro_block_size, ...
                                                     1:size(frame,2)/simulation_parameter.macro_block_size, ...
                                                     1:plane, ...
                                                     1:simulation_parameter.m);
        for k = 1:plane
            for i = 1:size(trancate_translation_temp,1)
                for j = 1:size(trancate_translation_temp,2)
                    for z = 1:simulation_parameter.m
                        y.element_plane_translation{i,j,k}(z,1) = trancate_translation_temp(i,j,k,z);
                    end
                end
            end
        end
        
        translation_temp_q = cell2mat(quantization_coded_sub_frame);
        %___Trancate to original size
        trancate_translation_temp_q = translation_temp_q(1:size(frame,1)/simulation_parameter.macro_block_size, ...
                                                     1:size(frame,2)/simulation_parameter.macro_block_size, ...
                                                     1:plane, ...
                                                     1:simulation_parameter.m);
        for k = 1:plane
            for i = 1:size(trancate_translation_temp_q,1)
                for j = 1:size(trancate_translation_temp_q,2)
                    for z = 1:simulation_parameter.m
                        y.element_plane_translation_q{i,j,k}(z,1) = trancate_translation_temp_q(i,j,k,z);
                    end
                end
            end
        end
        
        disp('Done');
        
        disp('Recovering...');
        for k = 1:plane
            for i = 1:size(frame,1)/simulation_parameter.macro_block_size
                for j = 1:size(frame,2)/simulation_parameter.macro_block_size
                     reconstructed_image{i,j,k} = BCS_reconstruction(y.element_plane_translation{i,j,k}, ...
                                                                     simulation_parameter.phi, ...
                                                                     simulation_parameter.theta, ...
                                                                     simulation_parameter.reconstruction_algorithm, ...
                                                                     simulation_parameter.macro_block_size);
                     res_reconstructed_image{i,j,k} = BCS_reconstruction(y.element_plane_translation_q{i,j,k}, ...
                                                                         simulation_parameter.phi, ...
                                                                         simulation_parameter.theta, ...
                                                                         simulation_parameter.reconstruction_algorithm, ...
                                                                         simulation_parameter.macro_block_size);
                end
            end
        end
        disp('Recovering Done');
        video{matrix_depth,frame_number}            = uint8(cell2mat(reconstructed_image));
        res_video{matrix_depth,frame_number}        = uint8(cell2mat(res_reconstructed_image));
        for k = 1:plane
            video{matrix_depth,frame_number}(:,:,k) = medfilt2(video{matrix_depth,frame_number}(:,:,k),[3 3], 'indexed');
            video{matrix_depth,frame_number}(:,:,k) = wdenoise2(video{matrix_depth,frame_number}(:,:,k),'CycleSpinning',1);
            
        end
        %___QUATITATIVE MATRICES___
        quantitative_metrices.psnr(frame_number)                     = psnr(video{matrix_depth,frame_number}, uint8(frame));
        quantitative_metrices.ssim(frame_number)                     = ssim(video{matrix_depth,frame_number}, uint8(frame));
        quantitative_metrices.rmse(frame_number)                     = sqrt(mean(((video{matrix_depth,frame_number}(:))-uint8(frame(:))).^2));
        quantitative_metrices.entropy_bpp_compressed(frame_number)   = entropy_bpp_compressed;
        imwrite(video{matrix_depth,frame_number}, strcat('final', num2str(z), '.png'))
    end
end

% profile report
% profile off