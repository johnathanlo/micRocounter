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
List findMS(std::string fileloc, NumericVector minrepeats, NumericVector tolerancefactors)
{
    char ch, bpmem[12], twomer[3], threemer[4], fourmer[5], fivemer[6], sixmer[7],
        header_str[200];
    const char *filename = fileloc.c_str();
    int i, j, rep_counter, skips = 0, flag = 0, oncount = 0, tot_microsat_content = 0,
        header = 0, header_length = 0, hloc = 0, trun = 0, gloc = 0, mbp_gloc = 0;
    struct merloc twom = {.length = (int *)malloc(sizeof(int)*1), .loc = (int *)malloc(sizeof(int)*1),
                          .tolerance = (int)tolerancefactors[0], .minrepeats = (int)minrepeats[0],
                          .seq = (char **)malloc(sizeof(char*)*1), .header_seqs = (char **)malloc(sizeof(char*)*1)},
                  threem = {.length = (int *)malloc(sizeof(int)*1), .loc = (int *)malloc(sizeof(int)*1),
                            .tolerance = (int)tolerancefactors[1], .minrepeats = (int)minrepeats[1],
                            .seq = (char **)malloc(sizeof(char*)*1), .header_seqs = (char **)malloc(sizeof(char*)*1)},
                  fourm = {.length = (int *)malloc(sizeof(int)*1), .loc = (int *)malloc(sizeof(int)*1),
                           .tolerance = (int)tolerancefactors[2], .minrepeats = (int)minrepeats[2],
                           .seq = (char **)malloc(sizeof(char*)*1), .header_seqs = (char **)malloc(sizeof(char*)*1)},
                  fivem = {.length = (int *)malloc(sizeof(int)*1), .loc = (int *)malloc(sizeof(int)*1),
                           .tolerance = (int)tolerancefactors[3], .minrepeats = (int)minrepeats[3],
                           .seq = (char **)malloc(sizeof(char*)*1), .header_seqs = (char **)malloc(sizeof(char*)*1)},
                  sixm= {.length = (int *)malloc(sizeof(int)*1), .loc = (int *)malloc(sizeof(int)*1),
                         .tolerance = (int)tolerancefactors[4], .minrepeats = (int)minrepeats[4],
                         .seq = (char **)malloc(sizeof(char*)*1), .header_seqs = (char **)malloc(sizeof(char*)*1)};
    clock_t begin = clock();
    FILE *gdat;
    gdat = fopen(filename, "r");

    while((ch = tolower(getc(gdat))) != EOF)
    {
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
      if((ch == 'a' || ch == 'g' || ch == 'c' || ch == 't') && header == 0)
      {
        for( i = 11; i>0; i --)
        {
          bpmem[i] = bpmem[i-1];
        }

        bpmem[0] = ch;

//if repeats are detected, these nested if statements will execute - in the future change these to a recursive function
        if(flag == 1)
        {
            oncount ++;
            if(twom.flag == 0)
            {
              if(threem.flag == 0)
              {
                if(fourm.flag == 0)
                {
                  if(fivem.flag == 0)
                  {
                    if(sixm.flag == 0)
                    {
                    }
                    else
                    {
                      if(oncount%6 == 0)
                      {
                        if(sixmer[0] == bpmem[6] && sixmer[1] == bpmem[7] && sixmer[2] == bpmem[8] && sixmer[3] == bpmem[9] && sixmer[4] == bpmem[10] && sixmer[5] == bpmem[11])
                        {
                          if(rep_counter == sixm.minrepeats)
                          {
                            sixm.totlocs++;//total loci of qualifying sixmers
                            sixm.length = (int*) realloc(sixm.length, sixm.totlocs * sizeof(int));//corresponding number of repeats for sixmer at that locus
                            sixm.length[sixm.totlocs-1] = sixm.minrepeats;//starting length of a qualifying sixmer = minimum length to qualify
                            sixm.loc = (int*) realloc(sixm.loc, sixm.totlocs * sizeof(int));
                            sixm.loc[sixm.totlocs-1] = hloc - 6*(sixm.minrepeats+1);//assign qualifying sixmer a starting location
                            sixm.seq = (char**) realloc(sixm.seq, sixm.totlocs* sizeof(char*));
                            sixm.seq[sixm.totlocs-1] = (char*) malloc(7*sizeof(char));
                            sixm.header_seqs = (char**) realloc(sixm.header_seqs, sixm.totlocs* sizeof(char*));
                            sixm.header_seqs[sixm.totlocs-1] = (char*) malloc(151*sizeof(char));
                            tot_microsat_content += 6*sixm.minrepeats;

                            //reverse the string
                            for(i=0, j = 5; j>=0; i++, j--)
                            {
                              sixm.seq[sixm.totlocs-1][i] = sixmer[j];
                            }
                            sixm.seq[sixm.totlocs-1][6] = '\0';

                            i = 0;
                            while(i < 149 && header_str[i] != '\0')
                            {
                              sixm.header_seqs[sixm.totlocs-1][i] = header_str[i];
                              i++;
                            }
                            sixm.header_seqs[sixm.totlocs-1][i] = '\0';

                          }

                          if(rep_counter > sixm.minrepeats)
                          {
                            sixm.length[sixm.totlocs-1] ++;
                            tot_microsat_content += 6;
                          }
                          rep_counter++;
                        }
                        else
                        {
                          skips++;
                          oncount = 5;
                          if(skips > sixm.tolerance)
                          {
                            flag = 0;
                            sixm.flag = 0;
                          }
                        }
                      }
                    }
                  }
                  else
                  {
                    if(oncount%5 == 0)
                    {
                      if(fivemer[0] == bpmem[5] && fivemer[1] == bpmem[6] && fivemer[2] == bpmem[7] && fivemer[3] == bpmem[8] && fivemer[4] == bpmem[9])
                      {
                        if(rep_counter == fivem.minrepeats)
                        {
                          fivem.totlocs++;//total loci of qualifying fivemers
                          fivem.length = (int*) realloc(fivem.length, fivem.totlocs * sizeof(int));//corresponding number of repeats for fivemer at that locus
                          fivem.length[fivem.totlocs-1] = fivem.minrepeats;//starting length of a qualifying fivemer = minimum length to qualify
                          fivem.loc = (int*) realloc(fivem.loc, fivem.totlocs * sizeof(int));
                          fivem.loc[fivem.totlocs-1] = hloc - 5*(fivem.minrepeats+1);//assign qualifying fivemer a starting location
                          fivem.seq = (char**) realloc(fivem.seq, fivem.totlocs* sizeof(char*));
                          fivem.seq[fivem.totlocs-1] = (char*) malloc(6*sizeof(char));
                          fivem.header_seqs = (char**) realloc(fivem.header_seqs, fivem.totlocs* sizeof(char*));
                          fivem.header_seqs[fivem.totlocs-1] = (char*) malloc(151*sizeof(char));
                          tot_microsat_content += 3*fivem.minrepeats;

                          //reverse the string
                          for(i=0, j = 4; j>=0; i++, j--)
                          {
                            fivem.seq[fivem.totlocs-1][i] = fivemer[j];
                          }
                          fivem.seq[fivem.totlocs-1][5] = '\0';

                          i = 0;
                          while(i < 149 && header_str[i] != '\0')
                          {
                            fivem.header_seqs[fivem.totlocs-1][i] = header_str[i];
                            i++;
                          }
                          fivem.header_seqs[fivem.totlocs-1][i] = '\0';
                        }

                        if(rep_counter > fivem.minrepeats)
                        {
                          fivem.length[fivem.totlocs-1] ++;
                          tot_microsat_content += 5;
                        }
                        rep_counter++;
                      }
                      else
                      {
                        skips ++;
                        oncount = 4;
                        if(skips>fivem.tolerance)
                        {
                          flag = 0;
                          fivem.flag = 0;
                        }
                      }
                    }
                  }
                }
                else
                {
                  if(oncount%4 == 0)
                  {
                    if(fourmer[0] == bpmem[4] && fourmer[1] == bpmem[5] && fourmer[2] == bpmem[6] && fourmer[3] == bpmem[7])
                    {
                      if(rep_counter == fourm.minrepeats)
                      {
                        fourm.totlocs++;//total loci of qualifying fourmers
                        fourm.length = (int*) realloc(fourm.length, fourm.totlocs * sizeof(int));//corresponding number of repeats for fourmer at that locus
                        fourm.length[fourm.totlocs-1] = fourm.minrepeats;//starting length of a qualifying fourmer = minimum length to qualify
                        fourm.loc = (int*) realloc(fourm.loc, fourm.totlocs * sizeof(int));
                        fourm.loc[fourm.totlocs-1] = hloc - 4*(fourm.minrepeats+1);//assign qualifying fourmer a starting location
                        fourm.seq = (char**) realloc(fourm.seq, fourm.totlocs* sizeof(char*));
                        fourm.seq[fourm.totlocs-1] = (char*) malloc(5*sizeof(char));
                        fourm.header_seqs = (char**) realloc(fourm.header_seqs, fourm.totlocs* sizeof(char*));
                        fourm.header_seqs[fourm.totlocs-1] = (char*) malloc(151*sizeof(char));
                        tot_microsat_content += 4*fourm.minrepeats;

                        //reverse the string
                        for(i=0, j = 3; j>=0; i++, j--)
                        {
                          fourm.seq[fourm.totlocs-1][i] = fourmer[j];
                        }
                        fourm.seq[fourm.totlocs-1][4] = '\0';

                        i = 0;
                        while(i < 149 && header_str[i] != '\0')
                        {
                          fourm.header_seqs[fourm.totlocs-1][i] = header_str[i];
                          i++;
                        }
                        fourm.header_seqs[fourm.totlocs-1][i] = '\0';
                      }

                      if(rep_counter > fourm.minrepeats)
                      {
                        fourm.length[fourm.totlocs-1] ++;
                        tot_microsat_content += 4;
                      }
                      rep_counter++;
                    }
                    else
                    {
                      oncount = 3; // reset oncount so that the next sequence automatically goes into this loop rather than skipping another x number of bases before going in again
                      skips++;
                      if(skips > fourm.tolerance)
                      {
                        flag = 0;
                        fourm.flag = 0;
                      }
                    }
                  }
                }
              }
              else
              {
                if(oncount%3 == 0)
                {
                  if(threemer[0] == bpmem[3] && threemer[1] == bpmem[4] && threemer[2] == bpmem[5])
                  {
                    if(rep_counter == threem.minrepeats)
                    {
                      threem.totlocs++;//total loci of qualifying threemers
                      threem.length = (int*) realloc(threem.length, threem.totlocs * sizeof(int));//corresponding number of repeats for threemer at that locus
                      threem.length[threem.totlocs-1] = threem.minrepeats;//starting length of a qualifying threemer = minimum length to qualify
                      threem.loc = (int*) realloc(threem.loc, threem.totlocs * sizeof(int));
                      threem.loc[threem.totlocs-1] = hloc - 3*(threem.minrepeats+1);//assign qualifying threemer a starting location
                      threem.seq = (char**) realloc(threem.seq, threem.totlocs* sizeof(char*));
                      threem.seq[threem.totlocs-1] = (char*) malloc(4*sizeof(char));
                      threem.header_seqs = (char**) realloc(threem.header_seqs, threem.totlocs* sizeof(char*));
                      threem.header_seqs[threem.totlocs-1] = (char*) malloc(151*sizeof(char));
                      tot_microsat_content += 3*threem.minrepeats;

                      //reverse the string
                      for(i=0, j = 2; j>=0; i++, j--)
                      {
                        threem.seq[threem.totlocs-1][i] = threemer[j];
                      }
                      threem.seq[threem.totlocs-1][3] = '\0';

                      i = 0;
                      while(i < 149 && header_str[i] != '\0')
                      {
                        threem.header_seqs[threem.totlocs-1][i] = header_str[i];
                        i++;
                      }
                      threem.header_seqs[threem.totlocs-1][i] = '\0';
                    }

                    if(rep_counter > threem.minrepeats)
                    {
                      threem.length[threem.totlocs-1] ++;
                      tot_microsat_content += 3;
                    }
                    rep_counter++;
                  }
                  else
                  {
                    skips++;
                    oncount = 2;
                    if(skips > threem.tolerance)
                    {
                      flag = 0;
                      threem.flag = 0;
                    }
                  }
                }
              }
            }
            else
            {
              if(oncount%2 == 0)
              {
                if(twomer[0] == bpmem[2] && twomer[1] == bpmem[3])
                {
                  if(rep_counter == twom.minrepeats)
                  {
                    twom.totlocs++;//total loci of qualifying twomers
                    twom.length = (int*) realloc(twom.length, twom.totlocs * sizeof(int));//corresponding number of repeats for twomer at that locus
                    twom.length[twom.totlocs-1] = twom.minrepeats;//starting length of a qualifying twomer = minimum length to qualify
                    twom.loc = (int*) realloc(twom.loc, twom.totlocs * sizeof(int));
                    twom.loc[twom.totlocs-1] = hloc - 2*(twom.minrepeats+1);//assign qualifying twomer a starting location
                    twom.seq = (char**) realloc(twom.seq, twom.totlocs* sizeof(char*));
                    twom.seq[twom.totlocs-1] = (char*) malloc(3*sizeof(char));
                    twom.header_seqs = (char**) realloc(twom.header_seqs, twom.totlocs* sizeof(char*));
                    twom.header_seqs[twom.totlocs-1] = (char*) malloc(151*sizeof(char));
                    tot_microsat_content += 2*twom.minrepeats;

                    //reverse the string
                    for(i=0, j = 1; j>=0; i++, j--)
                    {
                      twom.seq[twom.totlocs-1][i] = twomer[j];
                    }
                    twom.seq[twom.totlocs-1][2] = '\0';

                    i = 0;
                    while(i < 149 && header_str[i] != '\0')
                    {
                      twom.header_seqs[twom.totlocs-1][i] = header_str[i];
                      i++;
                    }
                    twom.header_seqs[twom.totlocs-1][i] = '\0';

                  }

                  if(rep_counter > twom.minrepeats)
                  {
                    twom.length[twom.totlocs-1] ++;
                    tot_microsat_content += 2;
                  }
                  rep_counter++;
                }
                else
                {
                  skips++;
                  oncount = 1;
                  if(skips >  twom.tolerance)
                  {
                    flag = 0;
                    twom.flag = 0;
                  }
                }
              }
            }
          }

        if(flag == 0)
        {
          if(bpmem[0] == bpmem[2] && bpmem[1] == bpmem[3])
          {
            flag = 1;
            twom.flag = 1;
            rep_counter = 2;
            oncount = 0;
            skips = 0;

            for(i=0; i<2; i++)
              twomer[i] = bpmem[i];
            twomer[2] = '\0';

          }
          else
          {
            if(bpmem[0] == bpmem[3] && bpmem[1] == bpmem[4] && bpmem[2] == bpmem[5])
            {
              flag = 1;
              threem.flag = 1;
              rep_counter = 2;
              oncount = 0;
              skips = 0;

              for(i=0; i<3; i++)
                threemer[i] = bpmem[i];
              threemer[3] = '\0';
            }
            else
            {
              if(bpmem[0] == bpmem[4] && bpmem[1] == bpmem[5] && bpmem[2] == bpmem[6] && bpmem[3] == bpmem[7])
              {
                flag = 1;
                fourm.flag = 1;
                rep_counter = 2;
                oncount = 0;
                skips = 0;

                for(i=0; i<4; i++)
                  fourmer[i] = bpmem[i];
                fourmer[4] = '\0';
              }
              else
              {
                if(bpmem[0] == bpmem[5] && bpmem[1] == bpmem[6] && bpmem[2] == bpmem[7] && bpmem[3] == bpmem[8] && bpmem[4] == bpmem[9])
                {
                  flag = 1;
                  fivem.flag = 1;
                  rep_counter = 2;
                  oncount = 0;
                  skips = 0;

                  for(i=0; i<5; i++)
                    fivemer[i] = bpmem[i];
                  fivemer[5] = '\0';
                }
                else
                {
                  if(bpmem[0] == bpmem[6] && bpmem[1] == bpmem[7] && bpmem[2] == bpmem[8] && bpmem[3] == bpmem[9] && bpmem[4] == bpmem[10] && bpmem[5] == bpmem[11])
                  {
                    flag = 1;
                    sixm.flag = 1;
                    rep_counter = 2;
                    oncount = 0;
                    skips = 0;

                    for(i=0; i<6; i++)
                      sixmer[i] = bpmem[i];
                    sixmer[6] = '\0';
                  }
                }
              }
            }
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

    mbp_gloc = mbp_gloc;


    NumericVector Rtwomlocs(twom.totlocs),
                  Rtwomlens(twom.totlocs),
                  Rthreemlocs(threem.totlocs),
                  Rthreemlens(threem.totlocs),
                  Rfourmlocs(fourm.totlocs),
                  Rfourmlens(fourm.totlocs),
                  Rfivemlocs(fivem.totlocs),
                  Rfivemlens(fivem.totlocs),
                  Rsixmlocs(sixm.totlocs),
                  Rsixmlens(sixm.totlocs);
    StringVector  Rtwomseqs(twom.totlocs),
                  Rtwomheaders(twom.totlocs),
                  Rthreemseqs(threem.totlocs),
                  Rthreemheaders(threem.totlocs),
                  Rfourmseqs(fourm.totlocs),
                  Rfourmheaders(fourm.totlocs),
                  Rfivemseqs(fivem.totlocs),
                  Rfivemheaders(fivem.totlocs),
                  Rsixmseqs(sixm.totlocs),
                  Rsixmheaders(sixm.totlocs);


    for(i=0; i<twom.totlocs; i++)
    {
      Rtwomlocs[i] = twom.loc[i];
      Rtwomlens[i] = twom.length[i];
      Rtwomseqs[i] = twom.seq[i];
      Rtwomheaders[i] = twom.header_seqs[i];
    }
    for(i=0; i<threem.totlocs; i++)
    {
      Rthreemlocs[i] = threem.loc[i];
      Rthreemlens[i] = threem.length[i];
      Rthreemseqs[i] = threem.seq[i];
      Rthreemheaders[i] = threem.header_seqs[i];
    }
    for(i=0; i<fourm.totlocs; i++)
    {
      Rfourmlocs[i] = fourm.loc[i];
      Rfourmlens[i] = fourm.length[i];
      Rfourmseqs[i] = fourm.seq[i];
      Rfourmheaders[i] = fourm.header_seqs[i];
    }
    for(i=0; i<fivem.totlocs; i++)
    {
      Rfivemlocs[i] = fivem.loc[i];
      Rfivemlens[i] = fivem.length[i];
      Rfivemseqs[i] = fivem.seq[i];
      Rfivemheaders[i] = fivem.header_seqs[i];
    }
    for(i=0; i<sixm.totlocs; i++)
    {
      Rsixmlocs[i] = sixm.loc[i];
      Rsixmlens[i] = sixm.length[i];
      Rsixmseqs[i] = sixm.seq[i];
      Rsixmheaders[i] = sixm.header_seqs[i];
    }

    List          Rtwom(List::create(Named("Loci") = Rtwomlocs, _["Lengths"] = Rtwomlens, _["Sequence"] = Rtwomseqs, _["SequenceNames"] = Rtwomheaders)),
                  Rthreem(List::create(Named("Loci") = Rthreemlocs, _["Lengths"] = Rthreemlens, _["Sequence"] = Rthreemseqs, _["SequenceNames"] = Rthreemheaders)),
                  Rfourm(List::create(Named("Loci") = Rfourmlocs, _["Lengths"] = Rfourmlens, _["Sequence"] = Rfourmseqs, _["SequenceNames"] = Rfourmheaders)),
                  Rfivem(List::create(Named("Loci") = Rfivemlocs, _["Lengths"] = Rfivemlens, _["Sequence"] = Rfivemseqs, _["SequenceNames"] = Rfivemheaders)),
                  Rsixm(List::create(Named("Loci") = Rsixmlocs, _["Lengths"] = Rsixmlens, _["Sequence"] = Rsixmseqs, _["SequenceNames"] = Rsixmheaders));

    return Rcpp::List::create(Rcpp::Named("Twomers") = Rtwom,
                              Rcpp::Named("Threemers") = Rthreem,
                              Rcpp::Named("Fourmers") = Rfourm,
                              Rcpp::Named("Fivemers") = Rfivem,
                              Rcpp::Named("Sixmers") = Rsixm,
                              Rcpp::Named("Genome Size (Mbp)") = mbp_gloc,
                              Rcpp::Named("Total Microsat Content") = tot_microsat_content);
    return 0;
}
