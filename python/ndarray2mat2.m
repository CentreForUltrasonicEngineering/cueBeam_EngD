% ndarray2mat2
function out=ndarray2mat2(in)
sizes=cell(in.shape); size1=sizes{1}.double; size2=sizes{2}.double;
out=reshape(typecast(uint8(in.base),'single'),size2*2,size1);
out=out(1:2:end-1,:)+1i*out(2:2:end,:);
out=out';
