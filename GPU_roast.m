subs = [9051 9021 9023 9031 9032 9047 9054]; %1892 101549 300100 300700 301051]; % 9054 300700
for sub = subs
   tic
   parfor e = 1:2485 % HARDCODED
       parallel.gpu.enableCUDAForwardCompatibility(true)
%        gpuDevice( 1 + mod( labindex - 1, gpuDeviceCount ) );
       try exhaustiveSearch_sGauss_GPU_III(sub,e)
       catch ME
           fprintf('sub-%04.f %01.f FAILED: %s\n',sub,e,ME.message)
       end
   end
   fprintf('sub-%04.f compelete: ~%.02f sec\n',sub,toc)
end
