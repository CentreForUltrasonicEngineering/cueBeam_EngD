% ndarray2mat2r
% real-valued version of the ndarray2mat2
function out=ndarray2mat2r(in)
sizes=cell(in.shape); size1=sizes{1}.double; size2=sizes{2}.double;
out=reshape(typecast(uint8(in.base),'single'),size2,size1);
%out=out(1:2:end-1,:)+1i*out(2:2:end,:);
out=out';
