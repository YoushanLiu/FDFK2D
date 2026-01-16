!=================================================================!
subroutine allocate_unit(iunit)

implicit none

integer, intent(inout) :: iunit


logical unit_opened


iunit = 21
inquire(unit = iunit, opened = unit_opened)
do while(unit_opened)
   iunit = iunit + 1
   inquire(unit = iunit, opened = unit_opened)
end do

end subroutine allocate_unit
!=================================================================!
subroutine elapsed_time(dtm_start, dtm_end, dtm_incr)

implicit none

integer(4), intent(inout) :: dtm_start(1:8)
integer(4), intent(inout) :: dtm_end(1:8)
integer(4), intent(inout) :: dtm_incr(1:8)


logical, external :: is_leap_year

integer(4) dtm_swp(1:8)

integer(2) months_start(1:12)
integer(2) months_end(1:12)

integer(2) iyear, year_min
integer(2) year_start, year_end, year_incr
integer(2) month_start, month_end, month_incr
integer(2) day_start, day_end, day_incr
integer(2) hour_start, hour_end, hour_incr
integer(2) min_start, min_end, min_incr
integer(2) sec_start, sec_end, sec_incr
integer(4) millisec_start, millisec_end, millisec_incr

integer(8) days, days_start, days_end

logical is_swap


dtm_start(6) = dtm_start(6) - dtm_start(4)
dtm_end(6) = dtm_end(6) - dtm_end(4)

if (dtm_start(1) > dtm_end(1)) then
   dtm_swp = dtm_start
   dtm_start = dtm_end
   dtm_end = dtm_swp
   is_swap = .true.
else
   is_swap = .false.
end if

if (.not.(is_swap) .and. (dtm_start(1) == dtm_end(1)) .and. (dtm_start(2) > dtm_end(2))) then
   dtm_swp = dtm_start
   dtm_start = dtm_end
   dtm_end = dtm_swp
   is_swap = .true.
else
   is_swap = .false.
end if

if (.not.(is_swap) .and. (dtm_start(2) == dtm_end(2)) .and. (dtm_start(3) > dtm_end(3))) then
   dtm_swp = dtm_start
   dtm_start = dtm_end
   dtm_end = dtm_swp
   is_swap = .true.
else
   is_swap = .false.
end if

if (.not.(is_swap) .and. (dtm_start(3) == dtm_end(3)) .and. (dtm_start(5) > dtm_end(5))) then
   dtm_swp = dtm_start
   dtm_start = dtm_end
   dtm_end = dtm_swp
   is_swap = .true.
else
   is_swap = .false.
end if

if (.not.(is_swap) .and. (dtm_start(5) == dtm_end(5)) .and. (dtm_start(6) > dtm_end(6))) then
   dtm_swp = dtm_start
   dtm_start = dtm_end
   dtm_end = dtm_swp
   is_swap = .true.
else
   is_swap = .false.
end if

if (.not.(is_swap) .and. (dtm_start(6) == dtm_end(6)) .and. (dtm_start(7) > dtm_end(7))) then
   dtm_swp = dtm_start
   dtm_start = dtm_end
   dtm_end = dtm_swp
   is_swap = .true.
else
   is_swap = .false.
end if

if (.not.(is_swap) .and. (dtm_start(7) == dtm_end(7)) .and. (dtm_start(8) > dtm_end(8))) then
   dtm_swp = dtm_start
   dtm_start = dtm_end
   dtm_end = dtm_swp
   is_swap = .true.
else
   is_swap = .false.
end if

year_start = dtm_start(1)
month_start = dtm_start(2)
day_start = dtm_start(3)
hour_start = dtm_start(5)
min_start = dtm_start(6)
sec_start = dtm_start(7)
millisec_start = dtm_start(8)

year_end = dtm_end(1)
month_end = dtm_end(2)
day_end = dtm_end(3)
hour_end = dtm_end(5)
min_end = dtm_end(6)
sec_end = dtm_end(7)
millisec_end = dtm_end(8)

months_start(1)  = 31; months_start(2)  = 28; months_start(3)  = 31
months_start(4)  = 30; months_start(5)  = 31; months_start(6)  = 30
months_start(7)  = 31; months_start(8)  = 31; months_start(9)  = 30
months_start(10) = 31; months_start(11) = 30; months_start(12) = 31

months_end = months_start

if (is_leap_year(year_start)) months_start(2) = months_start(2) + 1
if (is_leap_year(year_end)) months_end(2) = months_end(2) + 1

year_min = min(year_start,year_end)

days_start = 0
do iyear = year_min, year_start-1, 1

   days_start = days_start + 365
   if (is_leap_year(iyear)) days_start = days_start + 1

end do
days_start = days_start + sum(months_start(1:month_start-1))
days_start = days_start + day_start

days_end = 0
do iyear = year_min, year_end-1, 1

   days_end = days_end + 365
   if (is_leap_year(iyear)) days_end = days_end + 1

end do
days_end = days_end + sum(months_end(1:month_end-1))
days_end = days_end + day_end

days = abs(days_end - days_start)
year_incr = int(days / 365)
month_incr = int((days - 365*year_incr ) / 30)
day_incr = days - year_incr*365 - month_incr*30

hour_incr = hour_end - hour_start
min_incr = min_end - min_start
sec_incr = sec_end - sec_start
millisec_incr = millisec_end - millisec_start

if (millisec_incr < 0) then
   sec_incr = sec_incr - 1
   millisec_incr = millisec_incr + 1000
end if

if (sec_incr < 0) then
   min_incr = min_incr - 1
   sec_incr = sec_incr + 60
end if

if (min_incr < 0 ) then
   hour_incr = hour_incr - 1
   min_incr = min_incr + 60
end if

if (hour_incr < 0 ) then
   day_incr = day_incr - 1
   hour_incr = hour_incr + 24
end if

! output
dtm_incr(1) = year_incr
dtm_incr(2) = month_incr
dtm_incr(3) = day_incr
dtm_incr(4) = 0
dtm_incr(5) = hour_incr
dtm_incr(6) = min_incr
dtm_incr(7) = sec_incr
dtm_incr(8) = millisec_incr
dtm_start(6) = dtm_start(6) + dtm_start(4)
dtm_end(6) = dtm_end(6) + dtm_end(4)


end subroutine elapsed_time
!=================================================================!
logical function is_leap_year(iyear)

implicit none

integer(2), intent(in) :: iyear

if (((0 /= mod(iyear,100)) .and. (0 == mod(iyear,4)) ) .or. (0 == mod(iyear,400))) then
   is_leap_year = .true.
else
   is_leap_year = .false.
end if

return

end function is_leap_year
!=================================================================!





