CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ovrwrt (line, iwhich)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Taken from the DAOPhot code by Stetson.
C     This is a cool routine for writing lots of info to the screen
C     and keeping everything on the same line.
C
      character*(*) line
      character*79 output
      integer len
      if (iwhich .eq. 1) then
         write (6,1) line
    1    format (a)
      else if (iwhich .eq. 2) then
         if (len(line) .lt. 79) then
            output = ' '
            output = line
            write (6,2) output, char(13), char(13)
            write (6,2) output, char(13), char(13)
            write (6,2) output, char(13), char(13)
    2       format (a, 2a1, $)
         else
            write (6,2) line, char(13), char(13)
         end if
      else if (iwhich .eq. 3) then
         write (6,3) line
    3    format (a)
      else
         write (6,4) line, char(13), char(13)
    4    format (/a, 2a1, $)
         write (6,2) line, char(13), char(13)
      end if
      return
      end