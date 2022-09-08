function v = zigzag(u)
% returns a vector v with the elements of a matrix u in zigzag order.
% Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use 
% Redistributed with modifications made by Ruiyuan Lin on 2 July, 2015.
[M,N] = size(u);  m=1; n=1; v=zeros(M*N,1); v(1)=u(m,n);  d='r';
for i=2:M*N
   switch d
     case 'u',  m=m-(m>1); n=n+(n<N); v(i) = u(m,n);  
                if n==N,  d='d';  elseif m==1, d='r'; end
     case 'l',  m=m+(m<M); n=n-(n>1); v(i) = u(m,n);  
                if m==M, d='r'; elseif n==1, d='d'; end  
     case 'd',  m=m+(m<M); v(i) = u(m,n);  
                if n==1,  d='u';  else  d='l';  end
     case 'r',  n=n+(n<N); v(i) = u(m,n);  
                if m==1,  d='l';  else d='u';  end
   end
end
