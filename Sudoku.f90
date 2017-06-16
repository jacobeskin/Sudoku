program Sudoku
  !compilation: make -f makesudoku

  !Brute force sudoku solver. Solves sudoku puzzles with Monte Carlo Simulated
  !Annealing method. 

  !The program is run the usual way, ./Sudoku.exe, and it asks for the name of
  !the text file that has the puzzle and the name of the file it will create for
  !the solution. The program prints out the status of the solver at regular
  !intervals, giving the value of the parameter "counter" that is being 
  !minimized, runtime up to that point and number of steps. When the solution
  !is found (counter reaches 0) the same information is printed out.

  !As the puzzle is being read in by the program, it counts all the free
  !"writable" slots for each subregion (the 3x3 grids) and assigns an array of
  !derived data types for each subregion individually, the length of the arrays
  !being the number of empty slots in that subregion. Then it assigns a random
  !integer between 1 and 9 at each empty slot, and counts the total number of
  !mistakes in the puzzle into a variable "counter_0" (mistakes being dublicate
  !numbers in rows/columns and in the subregions, e.g. if there are three 9's
  !in one row, that counts as 2 mistakes). The initial "temperature" is set at
  !1000 and the simulated annealing process first lowers it linearly to almost
  !zero in one million steps and then rising it up to 90% of the previous
  !starting "temperature". This is continued until the "counter" variable
  !reaches zero or it seems like it has hit a wall, being the same after
  !multiple millions of iterations. A change of state is done by first
  !generating a random integer from 1 to 9 in order to choose the subregion,
  !and after this a random integer between 1 and the number of empty initial
  !slots depending on the subregion. Then as the writable slot has been
  !randomly chosen, the number in it is changed also randomly. After this the
  !counter is calculated again and compared to "counter_0", and the new state
  !is accepted or rejected in the usual way for simulated annealing routines. 

  !Input is a sudoku puzzle that is in the format of one of the exmaple files.
  !Output is the solved puzzle and time taken. Warning, code is not pretty! And
  !not all sudokus can be solved with this program, some very hard ones may
  !present a problem..

  !Used modules
  use mtdefs
  use mtmod
  implicit none

  ! Defining own variable type that holds the information of position and
  ! value of the writable slot in the sudoku puzzle
  type :: slot
     integer :: posx, posy, val
  end type slot

  ! All kinds of variables 
  integer :: puzzle(9,9), sub1(3,3), sub2(3,3), sub3(3,3)
  integer :: sub4(3,3), sub5(3,3), sub6(3,3), sub7(3,3), sub8(3,3), sub9(3,3)
  integer :: i, j, k, l, N_siman, x, y, iter_siman, iter_mc
  integer :: counter, counter_0, old_val, new_val, e
  integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9
  double precision :: u, start, end, T, T_0, end_mid
  character (len=80) :: filename, solutionfile
  ! Each subregion has its own array of slot-type variables
  type(slot), allocatable :: empty1(:), empty2(:), empty3(:)
  type(slot), allocatable :: empty4(:), empty5(:), empty6(:), empty7(:)
  type(slot), allocatable :: empty8(:), empty9(:)

  ! Initialize RNG
  call sgrnd(getseed(info=1))
  ! Begin clock counter for timing the code
  call cpu_time(start)

  ! Read in the the puzzle file
  print *,
  print *, "Name of the puzzle file:"
  read(*,*) filename
  filename = trim(filename)
  open(unit=1, file=filename, status='unknown')
  do i = 1,9
     read(1,*) puzzle(i,1:9)
  end do
  close(1)
  print *, "Name of the solution file:"
  read(*,*) solutionfile
  print *, "Thanks!"
  print *,

  !----------------------------------------------------
  ! Replace the zeros with random integers from 1 to 9:
  !----------------------------------------------------

  ! Assign the 9 subregions and allocate space for the 
  ! arrays that hold the writable slots
  sub1 = puzzle(1:3,1:3)
  n1 = count(sub1==0)
  allocate(empty1(n1))
  sub2 = puzzle(1:3,4:6)
  n2 = count(sub2==0)
  allocate(empty2(n2))
  sub3 = puzzle(1:3,7:9)
  n3 = count(sub3==0)
  allocate(empty3(n3))
  sub4 = puzzle(4:6,1:3)
  n4 = count(sub4==0)
  allocate(empty4(n4))
  sub5 = puzzle(4:6,4:6)
  n5 = count(sub5==0)
  allocate(empty5(n5))
  sub6 = puzzle(4:6,7:9)
  n6 = count(sub6==0)
  allocate(empty6(n6))
  sub7 = puzzle(7:9,1:3)
  n7 = count(sub7==0)
  allocate(empty7(n7))
  sub8 = puzzle(7:9,4:6)
  n8 = count(sub8==0)
  allocate(empty8(n8))
  sub9 = puzzle(7:9,7:9)
  n9 = count(sub9==0)
  allocate(empty9(n9))

  ! Assign a random integer between 1-9 into each
  ! writable slot one subreagion at a time 
  ! (variable 'i' goes through each subregion)
  do i=1,9
     if (i==1) then ! subregion 1
        k = 1
        do j=1,3
           do l=1,3
              if (sub1(l,j)==0) then
                 empty1(k)%posx = j
                 empty1(k)%posy = l
                 empty1(k)%val = igrnd(1,9)
                 sub1(l,j) = empty1(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==2) then !subregion 2
        k = 1
        do j=1,3
           do l=1,3
              if (sub2(l,j)==0) then
                 empty2(k)%posx = j
                 empty2(k)%posy = l
                 empty2(k)%val = igrnd(1,9)
                 sub2(l,j) = empty2(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==3) then ! you know the drill...
        k = 1
        do j=1,3
           do l=1,3
              if (sub3(l,j)==0) then
                 empty3(k)%posx = j
                 empty3(k)%posy = l
                 empty3(k)%val = igrnd(1,9)
                 sub3(l,j) = empty3(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==4) then
         k = 1
        do j=1,3
           do l=1,3
              if (sub4(l,j)==0) then
                 empty4(k)%posx = j
                 empty4(k)%posy = l
                 empty4(k)%val = igrnd(1,9)
                 sub4(l,j) = empty4(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==5) then
         k = 1
        do j=1,3
           do l=1,3
              if (sub5(l,j)==0) then
                 empty5(k)%posx = j
                 empty5(k)%posy = l
                 empty5(k)%val = igrnd(1,9)
                 sub5(l,j) = empty5(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==6) then
         k = 1
        do j=1,3
           do l=1,3
              if (sub6(l,j)==0) then
                 empty6(k)%posx = j
                 empty6(k)%posy = l
                 empty6(k)%val = igrnd(1,9)
                 sub6(l,j) = empty6(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==7) then
        k = 1 
        do j=1,3
           do l=1,3
              if (sub7(l,j)==0) then
                 empty7(k)%posx = j
                 empty7(k)%posy = l
                 empty7(k)%val = igrnd(1,9)
                 sub7(l,j) = empty7(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==8) then
         k = 1
        do j=1,3
           do l=1,3
              if (sub8(l,j)==0) then
                 empty8(k)%posx = j
                 empty8(k)%posy = l
                 empty8(k)%val = igrnd(1,9)
                 sub8(l,j) = empty8(k)%val
                 k = k+1
              end if
           end do
        end do
     else if (i==9) then
         k = 1
        do j=1,3
           do l=1,3
              if (sub9(l,j)==0) then
                 empty9(k)%posx = j
                 empty9(k)%posy = l
                 empty9(k)%val = igrnd(1,9)
                 sub9(l,j) = empty9(k)%val
                 k = k+1
              end if
           end do
        end do
     end if
  end do
        
              

  ! Status of the puzzel is determined by a variable named 'counter'. Goal is 
  ! to reach counter==0. Value of 'counter' is higher the more mistakes there
  ! are in the puzzle. First calculate the initial value of the counter called
  ! 'counter_0':
  
  counter_0 = 0
  ! Check all columns for duplicate numbers 
  do i=1,9
     do j=1,9
        counter_0 = counter_0+abs((count(puzzle(:,j)==i)-1))
     end do
  end do
  ! Check all rows for duplicate numbers
  do i=1,9
     do j=1,9
        counter_0 = counter_0+abs((count(puzzle(j,:)==i)-1))
     end do
  end do
  ! Check all subregions for duplicate numbers
  do i=1,9
     counter_0 = counter_0+abs((count(sub1==i)-1))
     counter_0 = counter_0+abs((count(sub2==i)-1))
     counter_0 = counter_0+abs((count(sub3==i)-1))
     counter_0 = counter_0+abs((count(sub4==i)-1))
     counter_0 = counter_0+abs((count(sub5==i)-1))
     counter_0 = counter_0+abs((count(sub6==i)-1))
     counter_0 = counter_0+abs((count(sub7==i)-1))
     counter_0 = counter_0+abs((count(sub8==i)-1))
     counter_0 = counter_0+abs((count(sub9==i)-1))
  end do

  !--------------------------------------
  ! Start the simulated annealing process
  !--------------------------------------

  T_0 = 5000         ! Initial "temperature"
  N_siman = 1000000  ! number of annealing steps
  iter_mc = 0        ! Counter for the number of annealing runs
  do
     iter_mc = iter_mc+1
     T = T_0
     iter_siman = 0  ! Counter for the number of annealing steps
     do i=1,N_siman
        iter_siman = iter_siman+1
        counter = 0
        ! Assign the 9 subregions
        sub1 = puzzle(1:3,1:3)
        sub2 = puzzle(1:3,4:6)
        sub3 = puzzle(1:3,7:9)
        sub4 = puzzle(4:6,1:3)
        sub5 = puzzle(4:6,4:6)
        sub6 = puzzle(4:6,7:9)
        sub7 = puzzle(7:9,1:3)
        sub8 = puzzle(7:9,4:6)
        sub9 = puzzle(7:9,7:9)

        ! Choose a subregion by random
        k = igrnd(1,9)
        
        ! If region 1
        if (k==1) then
          
           ! Change one number in this subregion
           e = igrnd(1,n1)           ! Choose writable slot
           x = empty1(e)%posx
           y = empty1(e)%posy
           old_val = empty1(e)%val   ! Save the old number
           new_val = igrnd(1,9)      ! Generate new number
           empty1(e)%val = new_val
           sub1(y,x) = empty1(e)%val ! Update subregion
           ! Update puzzle
           puzzle(1:3,1:3) = sub1

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregions for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
     
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 empty1(e)%val = old_val
                 sub1(y,x) = empty1(e)%val
                 puzzle(1:3,1:3) = sub1
              end if
           end if

        end if

        ! If region 2
        if (k==2) then 
           
           ! Change one number in this subregion
           e = igrnd(1,n2)
           x = empty2(e)%posx
           y = empty2(e)%posy
           old_val = empty2(e)%val
           new_val = igrnd(1,9)
           empty2(e)%val = new_val
           sub2(y,x) = empty2(e)%val
           ! Update puzzle
           puzzle(1:3,4:6) = sub2

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 ! If new state is not accepted then change back
                 empty2(e)%val = old_val
                 sub2(y,x) = empty2(e)%val 
                 puzzle(1:3,4:6) = sub2 
              end if
           end if
           
        end if

        
        ! If region 3
        if (k==3) then 
           
           ! Change one number in this subregion
           e = igrnd(1,n3)
           x = empty3(e)%posx
           y = empty3(e)%posy
           old_val = empty3(e)%val
           new_val = igrnd(1,9)
           empty3(e)%val = new_val
           sub3(y,x) = empty3(e)%val
           ! Update puzzle
           puzzle(1:3,7:9) = sub3
           
           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 empty3(e)%val = old_val
                 sub3(y,x) = empty3(e)%val
                 puzzle(1:3,7:9) = sub3
              end if
           end if
           
        end if

        
        ! If region 4
        if (k==4) then
           
           ! Change one number in this subregion
           e = igrnd(1,n4)
           x = empty4(e)%posx
           y = empty4(e)%posy
           old_val = empty4(e)%val
           new_val = igrnd(1,9)
           empty4(e)%val = new_val
           sub4(y,x) = empty4(e)%val
           ! Update the puzzle
           puzzle(4:6,1:3) = sub4

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 ! If new state is not accepted then change back
                 empty4(e)%val = old_val
                 sub4(y,x) = empty4(e)%val
                 puzzle(4:6,1:3) = sub4
              end if
           end if
           
        end if


        ! If region 5
        if (k==5) then 
           
           ! Change one number in this subregion
           e = igrnd(1,n5)
           x = empty5(e)%posx
           y = empty5(e)%posy
           old_val = empty5(e)%val
           new_val = igrnd(1,9)
           empty5(e)%val = new_val
           sub5(y,x) = empty5(e)%val
           ! Update puzzle
           puzzle(4:6,4:6) = sub5

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 ! If new state is not accepted then change back
                 empty5(e)%val = old_val
                 sub5(y,x) = empty5(e)%val
                 ! Update puzzle
                 puzzle(4:6,4:6) = sub5
              end if
           end if
           
        end if


        ! If region 6
        if (k==6) then 
           
           ! Change one number in this subregion
           e = igrnd(1,n6)
           x = empty6(e)%posx
           y = empty6(e)%posy
           old_val = empty6(e)%val
           new_val = igrnd(1,9)
           empty6(e)%val = new_val
           sub6(y,x) = empty6(e)%val
           ! Update puzzle
           puzzle(4:6,7:9) = sub6

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 ! If new state is not accepted then change back
                 empty6(e)%val = old_val
                 sub6(y,x) = empty6(e)%val
                 ! Update puzzle
                 puzzle(4:6,7:9) = sub6
              end if
           end if
           
        end if


        ! If region 7
        if (k==7) then 
           
           ! Change one number in this subregion
           e = igrnd(1,n7)
           x = empty7(e)%posx
           y = empty7(e)%posy
           old_val = empty7(e)%val
           new_val = igrnd(1,9)
           empty7(e)%val = new_val
           sub7(y,x) = empty7(e)%val
           ! Update puzzle
           puzzle(7:9,1:3) = sub7

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 ! If new state is not accepted then change back
                 empty7(e)%val = old_val
                 sub7(y,x) = empty7(e)%val
                 ! Update puzzle
                 puzzle(7:9,1:3) = sub7
              end if
           end if
           
        end if

        
        ! If region 8
        if (k==8) then
          
          ! Change one number in this subregion
          e = igrnd(1,n8)
          x = empty8(e)%posx
          y = empty8(e)%posy
          old_val = empty8(e)%val
          new_val = igrnd(1,9)
          empty8(e)%val = new_val
          sub8(y,x) = empty8(e)%val 
          ! Update puzzle
          puzzle(7:9,4:6) = sub8

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 ! If new state is not accepted then change back
                 empty8(e)%val = old_val
                 sub8(y,x) = empty8(e)%val 
                 ! Update puzzle
                 puzzle(7:9,4:6) = sub8
              end if
           end if
           
        end if


        ! If region 9
        if (k==9) then
           
           ! Change one number in this subregion
           e = igrnd(1,n9)
           x = empty9(e)%posx
           y = empty9(e)%posy
           old_val = empty9(e)%val
           new_val = igrnd(1,9)
           empty9(e)%val = new_val
           sub9(y,x) = empty9(e)%val 
           ! Update puzzle
           puzzle(7:9,7:9) = sub9

           ! Calculate new counter and compare it to counter_0
           ! Check all columns for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(:,j)==l)-1))
              end do
           end do
           ! Check all rows for duplicate numbers
           do l=1,9
              do j=1,9
                 counter = counter+abs((count(puzzle(j,:)==l)-1))
              end do
           end do
           ! Check subregion for duplicates
           do l=1,9
              counter = counter+abs((count(sub1==l)-1))
              counter = counter+abs((count(sub2==l)-1))
              counter = counter+abs((count(sub3==l)-1))
              counter = counter+abs((count(sub4==l)-1))
              counter = counter+abs((count(sub5==l)-1))
              counter = counter+abs((count(sub6==l)-1))
              counter = counter+abs((count(sub7==l)-1))
              counter = counter+abs((count(sub8==l)-1))
              counter = counter+abs((count(sub9==l)-1))
           end do
              
           ! Check if we keep the change
           if ((counter-counter_0)<=0) then
              counter_0 = counter
           else if ((counter-counter_0)>0) then
              u = grnd()
              if (u<exp((counter_0-counter)/T)) then
                 counter_0 = counter
              else
                 ! If new state is not accepted then change back
                 empty9(e)%val = old_val
                 sub9(y,x) = empty9(e)%val 
                 ! Update puzzle
                 puzzle(7:9,7:9) = sub9
              end if
           end if
           
        end if

        ! If solution reached
        if (counter_0==0) then
           exit
        end if
        
        ! Lower "temperature"
        T = T - 0.00000001*T_0

     end do

     ! If solution reached
     if (counter_0==0) then
        exit
     end if

     ! Print some intermediate results
     call cpu_time(end_mid)
     if (mod(iter_mc,10)==1) then
        print *,
        print *, "Counter is: ", counter_0
        print *, "and runtime is:", (end_mid-start)/60
        print *, "Number of steps:", iter_mc*N_siman
        print *,
     end if

     ! Raise "temperature"
     T_0 = 0.9*T_0

  end do

  ! Write solution into a file
  open(unit=2, file=solutionfile, status='unknown')
  do i=1,9
     write(2,*) puzzle(i,:)
     write(2,*)
  end do
  close(2)

  ! Final announcement
  print *,
  print *, "Solution achieved in ", (iter_mc*N_siman)+iter_siman, " steps!"
  print *, "Counter_0 ", counter_0
        
  call cpu_time(end)

  ! Total elapsed time
  print *, "Elapsed time: ", (end-start)/60
  print *, 

end program Sudoku

  
