% python ndarray to matlab
function out=ndarray2mat(in)
sizes=cell(in.shape); size1=sizes{1}.double; size2=sizes{2}.double;
in_list=cell(in.tolist); % 2D array
out=zeros(size1,size2);
for idx1=1:size1
    cq1=cell(in_list{idx1});
    out(idx1,:)=cell2mat(cq1);
end


