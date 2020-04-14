function [S] = EMB_pdf2entropy(pdf)
%EMB_PDF2ENTROPY computes Shannon entropy of a pdf
l=log2(pdf);
l(l==-inf)=0;
S=-sum(pdf.*l,'all');
end

