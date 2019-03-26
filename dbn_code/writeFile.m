%Author: Daniel Ruiz-Perez
%Advisor Profs: Giri Narasimhan and Ziv Bar-Joseph
%Description: Given a filename, a matrix and the number of rows and columns, it writes the matrix to a file

function[] = writeFile(fName, cols, numdata, rows)
        fid = fopen(fName,'wt');
         for a = (1:(length(cols)-1))   
             fprintf(fid,  cols{a});
             fprintf(fid, '\t');
         end
         fprintf(fid,  cols{length(cols)});

         fprintf(fid, '\n');
         for a = rows 
                fprintf(fid,'%.1f', numdata(a,1));
                fprintf(fid, '\t');
                fprintf(fid,'%d', numdata(a,2));
                fprintf(fid, '\t');
                fprintf(fid,'%d', numdata(a,3));
                fprintf(fid, '\t');
            for b = (4:(length(numdata(1,:))-1))  
                 fprintf(fid,'%.4f', numdata(a,b));
                 fprintf(fid, '\t');
            end 
                fprintf(fid,'%.4f', numdata(a,length(numdata(1,:))));
                fprintf(fid, '\n');
         end
        fclose(fid);
        fclose all;
end

