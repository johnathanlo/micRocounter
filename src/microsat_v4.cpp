#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;


//structure to record numbers and locations of microsats
struct merloc
{
    int *length;
    int *loc;
    int tolerance;
    int minrepeats;
    char **seq;
    char **header_seqs;
    int totlocs, flag;

};
//[[Rcpp::export]]
List FindMS2(std::string fileloc, int replen, int minrepeats, int tolfac)
{
  char ch, *bp1, *bp2, *xmer, header_str[200];
  const char *filename = fileloc.c_str();
  int i, j, k, rep_counter = 0, skips = 0, flag = 0, oncount = 0, tot_microsat_content = 0,
      header = 0, header_length = 0, hloc = 0, trun = 0, gloc = 0, mbp_gloc = 0, match = 0;
  struct merloc xm = {.length = (int *)malloc(sizeof(int)*1), .loc = (int *)malloc(sizeof(int)*1),
                        .tolerance = tolfac, .minrepeats = minrepeats,
                        .seq = (char **)malloc(sizeof(char*)*1), .header_seqs = (char **)malloc(sizeof(char*)*1)};
  bp1 = (char *)malloc(sizeof(char)*replen);
  bp2 = (char *)malloc(sizeof(char)*replen);
  xmer = (char *)malloc(sizeof(char)*replen);
  clock_t begin = clock();
  FILE *gdat;
  gdat = fopen(filename, "r");

  while((ch = tolower(getc(gdat))) != EOF)
  {
/////scaff header stuff
    if(gloc > 1000000)
    {
      gloc = 0;
      mbp_gloc ++;
    }
    if(ch == '>')
    {
      header = 1;
      header_length = 1;
      hloc = 0;
    }
    if(ch != '\n' && ch != '>' && header == 1)
    {
      header_str[header_length - 1] = ch;
      header_length ++;
    }
    if(ch == '\n' && header == 1)
    {
      header = 0;
      header_str[header_length] = '\0';
      if (header_length > 150)
      {
        trun = 1;
      }
    }
    if(ch != '\n' && header == 0)
    {
      gloc++;
      hloc ++;
    }
///end scaff header stuff

///start scanning genome file
    if((ch == 'a' || ch == 'g' || ch == 'c' || ch == 't') && header == 0)
    {
///update character buffers
      for( i = 0; i<replen-1; i ++)
      {
        bp2[i] = bp2[i+1];
      }
      bp2[replen-1] = bp1[0];
      for(i = 0; i<replen-1; i++)
      {
        bp1[i] = bp1[i+1];
      }
      bp1[replen-1] = ch;
///finish updating char buffers

///if already in a microsatellite, check that matches continue
      if(flag == 1)
      {
        oncount++;
        if(oncount%replen == 0)
        {
          match = 1;
          for(i= 0; i < replen; i++)
          {
            if(bp1[i] != xmer[i])
            {
              match = 0;
            }
          }
          if(match == 1)
          {
            rep_counter++;
            if(rep_counter == xm.minrepeats)
            {
              xm.totlocs++;
              xm.length = (int*) realloc(xm.length, xm.totlocs * sizeof(int));//corresponding number of repeats for sixmer at that locus
              xm.length[xm.totlocs-1] = xm.minrepeats;//starting length of a qualifying sixmer = minimum length to qualify
              xm.loc = (int*) realloc(xm.loc, xm.totlocs * sizeof(int));
              xm.loc[xm.totlocs-1] = hloc - replen*(xm.minrepeats);//assign qualifying sixmer a starting location
              xm.seq = (char**) realloc(xm.seq, xm.totlocs* sizeof(char*));
              xm.seq[xm.totlocs-1] = (char*) malloc((replen+1)*sizeof(char));
              xm.header_seqs = (char**) realloc(xm.header_seqs, xm.totlocs* sizeof(char*));
              xm.header_seqs[xm.totlocs-1] = (char*) malloc(151*sizeof(char));
              tot_microsat_content += replen*xm.minrepeats;

//copy string over
              for(i=0;i<replen;i++)
              {
                xm.seq[xm.totlocs-1][i] = xmer[i];
              }
              xm.seq[xm.totlocs-1][replen] = '\0';

              i = 0;
              while(i < 149 && header_str[i] != '\0')
              {
                xm.header_seqs[xm.totlocs-1][i] = header_str[i];
                i++;
              }
              xm.header_seqs[xm.totlocs-1][i] = '\0';
            }

            if(rep_counter > xm.minrepeats)
            {
              xm.length[xm.totlocs-1] ++;
              tot_microsat_content += replen;
            }
          }

          else
          {
            skips++;
            oncount = replen-1;
            if(skips > xm.tolerance)
            {
              flag = 0;
              xm.flag = 0;
            }
          }
        }
      }
///if not in a microsatellite, check for matches given updated sequence
      if(flag == 0)
      {
        match = 1;
        for(i = 0; i<replen; i++)
        {
          if(bp1[i]!=bp2[i])
          {
            match = 0;
          }
        }

        if(match == 1 && replen >3)
        {
          for(i = 2; i< ((replen+1)/2 + (replen+1)%2); i++)
          {
            if(replen%i == 0)
            {
              match = 0;
              for(j=0; j<((replen/i)-1); j++)
              {
                for(k=0; k<i; k++)
                {
                  if(bp1[k]!=bp1[i+i*j+k])
                  {
                    match = 1;
                  }
                }
              }
              if(match == 0)
              {
                i = (replen+1)/2;
              }
            }
          }
        }

        if(match == 1)
        {
          flag = 1;
          rep_counter = 2;
          oncount = 0;
          skips = 0;

          for(i=0; i<replen; i++)
            xmer[i] = bp1[i];

          xmer[replen] = '\0';
        }
      }
    }
  }

  fclose(gdat);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("\nFile: %s ,Time: %lf", filename, time_spent);

  if(trun == 1)
  {
    printf("\nSequence names truncated to 150 characters.");
  }

  NumericVector Rxmlocs(xm.totlocs),
                Rxmlens(xm.totlocs);
  StringVector  Rxmseqs(xm.totlocs),
                Rxmheaders(xm.totlocs);


  for(i=0; i<xm.totlocs; i++)
  {
    Rxmlocs[i] = xm.loc[i];
    Rxmlens[i] = xm.length[i];
    Rxmseqs[i] = xm.seq[i];
    Rxmheaders[i] = xm.header_seqs[i];
  }

  List          Rxm(List::create(Named("Loci") = Rxmlocs, _["Lengths"] = Rxmlens, _["Sequence"] = Rxmseqs, _["SequenceNames"] = Rxmheaders));

  return Rcpp::List::create(Rcpp::Named("Xmers") = Rxm,
                            Rcpp::Named("Genome Size (Mbp)") = mbp_gloc,
                            Rcpp::Named("Total Microsat Content") = tot_microsat_content);
  return 0;

}
