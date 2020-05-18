function [S] = EMB_pdf2entropy(pdf,dim)
%EMB_PDF2ENTROPY computes Shannon entropy of a pdf
l=log2(pdf);
l(l==-inf)=0;
S=-EMB_sum(pdf.*l,dim);
end

