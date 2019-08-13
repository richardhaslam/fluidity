!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General PublicS
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module particles

  use fldebug
  use iso_c_binding, only: C_NULL_CHAR
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, &
       & PYTHON_FUNC_LEN, integer_size, real_size
  use futils, only: int2str, free_unit
  use elements
  use mpi_interfaces
  use parallel_tools
  use spud
  use parallel_fields
  use fields
  use state_module
  use field_options
  use detector_data_types
  use pickers
  use detector_tools
  use detector_parallel
  use detector_move_lagrangian
  use diagnostic_variables, only: field_tag, initialise_constant_diagnostics
  implicit none

  private

  public :: initialise_particle_positions, move_particles, write_particles_loop, destroy_particles, &
            update_particle_attributes_and_fields, checkpoint_particles_loop

  type(detector_linked_list), allocatable, dimension(:), target, save :: particle_lists !!Particle lists with dimension equal to the number of particle subgroups

contains

  subroutine initialise_particle_positions(filename,state)
    !!Initialise particles and set up particle file headers (per particle array)
    character(len = *), intent(in) :: filename
    type(state_type), dimension(:), intent(in) :: state

    character(len=FIELD_NAME_LEN) :: subname
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    type(vector_field), pointer :: xfield
    real :: current_time

    integer, dimension(3) :: attribute_size
    integer :: sub_particles
    integer :: i, k
    integer :: dim, particle_groups, total_arrays, list_counter
    integer, dimension(:), allocatable :: particle_arrays
    integer :: totaldet_global
    integer :: n_oldfields, phase, f
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield

    logical :: from_file

    ewrite(2,*), "In initialise_particles"

    !Check whether there are any particles.
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    total_arrays = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
       total_arrays = total_arrays + particle_arrays(i)
    end do
    allocate(particle_lists(total_arrays))

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)

    list_counter=1

    !Number of old_fields stored on particles
    n_oldfields = 0
    do phase = 1,size(state)
       do f = 1, size(state(phase)%scalar_names)
          sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(f))
          if (sfield%option_path=="".or.aliased(sfield)) then
             cycle
          else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
             n_oldfields = n_oldfields+1
          end if
       end do
       do f = 1, size(state(phase)%vector_names)
          vfield => extract_vector_field(state(phase),state(phase)%vector_names(f))
          if (vfield%option_path=="".or.aliased(vfield)) then
             cycle
          else if (have_option(trim(complete_field_path(vfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
             n_oldfields = n_oldfields+dim
          end if
       end do
       if (associated(state(phase)%tensor_names)) then
          do f = 1, size(state(phase)%tensor_names)
             tfield => extract_tensor_field(state(phase),state(phase)%tensor_names(f))
             if (tfield%option_path=="".or.aliased(tfield)) then
                cycle
             else if (have_option(trim(complete_field_path(tfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
                n_oldfields = n_oldfields+(dim**2)
             end if
          end do
       end if
    end do
    
    
    do i = 1,particle_groups
       !Set group_path
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          !subgroup_path
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          ! If the option
          ! "from_file" exists, it means we are continuing the simulation
          ! after checkpointing and the reading of the particle positions must be
          ! done from a file

          from_file=have_option(trim(subgroup_path)//"/initial_position/from_file")
          call get_option(trim(subgroup_path)//"/number_of_particles", sub_particles) !number of particles in subgroup
          call get_option(trim(subgroup_path)//"/name", subname) !Name of particle subgroup
          particle_lists(list_counter)%total_num_det=sub_particles
          allocate(particle_lists(list_counter)%detector_names(sub_particles))
          !Register this I/O list with a global list of detectors/particles
          call register_detector_list(particle_lists(list_counter))
          
          !Find number of attributes 
          attribute_size(:)=0
          if (have_option(trim(subgroup_path) // "/attributes")) then
             attribute_size(1)=option_count(trim(subgroup_path) // "/attributes/scalar_attribute")
             attribute_size(1)=attribute_size(1)+dim*option_count(trim(subgroup_path) // "/attributes/vector_attribute")
             attribute_size(1)=attribute_size(1)+(dim**2)*option_count(trim(subgroup_path) // "/attributes/tensor_attribute")
          end if
          if (option_count(trim(subgroup_path) // "/attributes/scalar_attribute/python_fields")>0 .or. &
             option_count(trim(subgroup_path) // "/attributes/vector_attribute/python_fields")>0 .or. &
             option_count(trim(subgroup_path) // "/attributes/tensor_attribute/python_fields")>0) then
             attribute_size(3) = n_oldfields
             attribute_size(2) = option_count(trim(subgroup_path) //"/attributes/scalar_attribute/python_fields/store_old_attribute")
             attribute_size(2) = attribute_size(2)+dim*option_count(trim(subgroup_path) //"/attributes/vector_attribute/python_fields/store_old_attribute")
             attribute_size(2) = attribute_size(2)+(dim**2)*option_count(trim(subgroup_path) //"/attributes/tensor_attribute/python_fields/store_old_attribute")
          end if
          
          ! Enable particles to drift with the mesh
          if (have_option("/particles/move_with_mesh")) then
             particle_lists(list_counter)%move_with_mesh=.true.
          end if
          
          ! Set flag for NaN particle output
          if (have_option("/particles/write_nan_outside_domain")) then
             particle_lists(list_counter)%write_nan_outside=.true.
          end if
          
          !Read particles from options
          if (.not.from_file) then
             
             call read_particles_from_python(sub_particles, subname, current_time, state, attribute_size, xfield, dim, subgroup_path, particle_lists(list_counter))
             
          else !Read particles from file
             
             call read_particles_from_file(sub_particles, subname, attribute_size, xfield, dim, subgroup_path, particle_lists(list_counter))
             
          end if !Read particles

          !call set_particle_output_file(sub_particles, subname, filename, attribute_size, xfield, subgroup_path, particle_lists(list_counter))
          !Get options for lagrangian particle movement
          call read_detector_move_options(particle_lists(list_counter), "/particles")
          list_counter = list_counter + 1
       end do
    end do

    ! And finally some sanity checks
    list_counter=1
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1,particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          call get_option(trim(subgroup_path)//"/name",subname)
          totaldet_global=particle_lists(list_counter)%length
          call allsum(totaldet_global)
          ewrite(2,*) "Found", particle_lists(list_counter)%length, "local and ", totaldet_global, "global particles for particle array ", trim(subname)
          assert(totaldet_global==particle_lists(list_counter)%total_num_det)
          list_counter = list_counter + 1
       end do
    end do

    deallocate(particle_arrays)
    
  end subroutine initialise_particle_positions

  subroutine read_particles_from_python(sub_particles, subname, current_time, state, attribute_size, xfield, dim, subgroup_path, p_list)
    ! Reading particles from a python function

    type(state_type), dimension(:), intent(in) :: state
    type(vector_field), pointer, intent(in) :: xfield
    type(detector_linked_list), intent(inout) :: p_list
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    character(len=FIELD_NAME_LEN), intent(in) :: subname
    real, intent(in) :: current_time
    integer, intent(in) :: sub_particles, dim
    integer, dimension(3), intent(in) :: attribute_size

    character(len=PYTHON_FUNC_LEN) :: func
    character(len = FIELD_NAME_LEN) :: fmt
    character(len=FIELD_NAME_LEN) :: particle_name
    real, allocatable, dimension(:,:) :: coords !array to hold coordinates of particles for initialisation
    real :: dt
    integer :: l, str_size

    ewrite(2,*) "Reading particles from options"
    call get_option(trim(subgroup_path)//"/initial_position/python", func)
    call get_option("/timestepping/timestep", dt)
    allocate(coords(dim,sub_particles))
    call set_detector_coords_from_python(coords, sub_particles, func, current_time)

    str_size=len_trim(int2str(sub_particles))
    fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

    do l=1,sub_particles
       write(particle_name, fmt) trim(subname)//"_", l
       call create_single_particle(p_list, xfield, coords(:,l), &
            l, LAGRANGIAN_DETECTOR, trim(particle_name), attribute_size)
    end do
    if (attribute_size(1)/=0) then
       call update_particle_subgroup_attributes_and_fields(state, current_time, dt, subgroup_path, p_list)
    end if
    
    deallocate(coords)

  end subroutine read_particles_from_python

  subroutine read_particles_from_file(sub_particles, subname, attribute_size, xfield, dim, subgroup_path, p_list)
    ! Reading from a binary file where the user has placed the particle information

    ! If reading from file:
    ! Particles checkpoint file names end in _par, with.groups appended for the header file
    ! and .attributes.dat appended for the binary data file that holds the positions and attributes

    type(vector_field), pointer, intent(in) :: xfield
    type(detector_linked_list), intent(inout) :: p_list
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    character(len=FIELD_NAME_LEN), intent(in) :: subname
    integer, intent(in) :: sub_particles, dim
    integer, dimension(3), intent(in) :: attribute_size

    real, allocatable, dimension(:,:) :: attribute_vals !array to hold particle attribute values for initialisation
    real, allocatable, dimension(:) :: packed_buff !array to hold particle attributes if from checkpoint file
    real, allocatable, dimension(:) :: particle_location !array to hold particle coordinates if from checkpoint file

    character(len = OPTION_PATH_LEN) :: particles_cp_filename
    character(len = FIELD_NAME_LEN) :: fmt
    character(len=FIELD_NAME_LEN) :: particle_name
    integer :: particle_checkpoint_unit=0
    integer :: l, m, str_size

    ewrite(2,*) "Reading particles from file"

    allocate(particle_location(dim))
    allocate(attribute_vals(3,maxval(attribute_size)))
    
    particle_checkpoint_unit=free_unit()
    call get_option(trim(subgroup_path) // "/initial_position/from_file/file_name",particles_cp_filename)
    
#ifdef STREAM_IO
    open(unit = particle_checkpoint_unit, file = trim(particles_cp_filename) //'.attributes.dat', &
         & action = "read", access = "stream", form = "unformatted")
#else
    FLExit("No stream I/O support required for particle checkpoints")
#endif
    
    !Read in particle locations from file       
    allocate(packed_buff(dim+sum(attribute_size)))
    str_size=len_trim(int2str(sub_particles))
    fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"
    
    do m=1,sub_particles
       write(particle_name, fmt) trim(subname)//"_", m
       read(particle_checkpoint_unit) packed_buff
       particle_location=packed_buff(1:dim)
       if (attribute_size(1)/=0) then
          attribute_vals(1,1:attribute_size(1))=packed_buff(dim+1:dim+attribute_size(1))
          attribute_vals(2,1:attribute_size(2))=packed_buff(dim+1+attribute_size(1):dim+attribute_size(1)+attribute_size(2))
          attribute_vals(3,1:attribute_size(3))=packed_buff(dim+1+attribute_size(1)+attribute_size(2):dim+sum(attribute_size))
       end if
       call create_single_particle(p_list, xfield, &
            particle_location, m, LAGRANGIAN_DETECTOR, trim(particle_name),attribute_size,attribute_vals= attribute_vals)
    end do
    deallocate(packed_buff)
    deallocate(attribute_vals)
    ewrite(2,*) "Finished read_particles_from_checkpoint"
    
  end subroutine read_particles_from_file

  subroutine set_particle_output_file(sub_particles, subname, filename, attribute_size, xfield, subgroup_path, p_list)
    !Sets up the particle output file type and output file.

    type(vector_field), pointer, intent(in) :: xfield
    type(detector_linked_list), intent(inout) :: p_list
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    character(len = *), intent(in) :: filename
    integer, intent(in) :: sub_particles
    integer, dimension(3), intent(in) :: attribute_size
    character(len=FIELD_NAME_LEN), intent(in) :: subname

    character(len=FIELD_NAME_LEN) :: attname
    character(len = OPTION_PATH_LEN) :: buffer
    integer :: column, ierror, m, n
    
    p_list%binary_output = .true.
    if (have_option("/particles/ascii_output")) then
       p_list%binary_output= .false.
       if(isparallel()) then
          FLExit("Error: No support for ascii detector output in parallel. Please use binary output by turning off the ascii_output option.")
       end if
    end if
    !!!!!Currently disabled until particle IO is revamped.
    
!!$    
!!$          ! Only the first process should write the header file
!!$    if (getprocno() == 1) then
!!$       p_list%output_unit=free_unit()
!!$       open(unit=p_list%output_unit, file=trim(filename)//'.particles.'//trim(subname), action="write")
!!$       
!!$       write(p_list%output_unit, '(a)') "<header>"
!!$       call initialise_constant_diagnostics(p_list%output_unit, &
!!$            binary_format = p_list%binary_output)
!!$       
!!$       ! Initial columns are elapsed time and dt.
!!$       buffer=field_tag(name="ElapsedTime", column=1, statistic="value")
!!$       write(p_list%output_unit, '(a)') trim(buffer)
!!$       buffer=field_tag(name="dt", column=2, statistic="value")
!!$       write(p_list%output_unit, '(a)') trim(buffer)
!!$       
!!$       ! Next columns contain the positions of all the particles
!!$       column=2
!!$       positionloop: do m=1,sub_particles
!!$          buffer=field_tag(name=p_list%detector_names(m), column=column+1,&
!!$               statistic="position", components=xfield%dim)
!!$          write(p_list%output_unit, '(a)') trim(buffer)
!!$          column=column+xfield%dim 
!!$       end do positionloop
!!$       
!!$       ! Next columns contain attributes of particles
!!$       attributeloop: do m=1,sub_particles
!!$          do n=1,attribute_size(1)
!!$             call get_option(trim(subgroup_path) // "/attributes/attribute["//int2str(n-1)//"]/name", attname)
!!$             buffer=field_tag(name=p_list%detector_names(m), column=column+1,&
!!$                  statistic=attname, components=1)
!!$             write(p_list%output_unit, '(a)')trim(buffer)
!!$             column=column+1
!!$          end do
!!$       end do attributeloop
!!$    end if
!!$    
!!$    if (getprocno() == 1) then
!!$       write(p_list%output_unit, '(a)') "</header>"
!!$       flush(p_list%output_unit)
!!$       ! when using mpi_subroutines to write into the particles file we need to close the file since 
!!$       ! filename.particles.dat needs to be open now with MPI_OPEN
!!$       if (p_list%binary_output) then
!!$          close(p_list%output_unit)
!!$       end if
!!$    end if
!!$    
!!$    
!!$    if (p_list%binary_output) then
!!$       ! bit of hack to delete any existing .particles.dat file
!!$       ! if we don't delete the existing .particles.dat would simply be opened for random access and 
!!$       ! gradually overwritten, mixing particle output from the current with that of a previous run
!!$       call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.particles.'//trim(subname)//'.dat', MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, p_list%mpi_fh, IERROR)
!!$       call MPI_FILE_CLOSE(p_list%mpi_fh, IERROR)
!!$       call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.particles.'//trim(subname)//'.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, p_list%mpi_fh, IERROR)
!!$       assert(ierror == MPI_SUCCESS)
!!$    end if

  end subroutine set_particle_output_file

  subroutine create_single_particle(detector_list, xfield, position, id, type, name, attribute_size, attribute_vals)
    ! Allocate a single particle, populate and insert it into the given list
    ! In parallel, first check if the particle would be local and only allocate if it is
    type(detector_linked_list), intent(inout) :: detector_list
    type(vector_field), pointer, intent(in) :: xfield
    real, dimension(xfield%dim), intent(in) :: position
    integer, intent(in) :: id, type
    character(len=*), intent(in) :: name

    type(detector_type), pointer :: detector
    type(element_type), pointer :: shape
    real, dimension(xfield%dim+1) :: lcoords
    integer :: element
    real, dimension(:,:), intent(in), optional :: attribute_vals
    integer, dimension(3), intent(in) :: attribute_size

    real ::  dt

    shape=>ele_shape(xfield,1)
    assert(xfield%dim+1==local_coord_count(shape))
    detector_list%detector_names(id)=name
    ! Determine element and local_coords from position
    call picker_inquire(xfield,position,element,local_coord=lcoords,global=.true.)
    call get_option("/timestepping/timestep", dt)
    ! If we're in parallel and don't own the element, skip this particle
    if (isparallel()) then
       if (element<0) return
       if (.not.element_owned(xfield,element)) return
    else
       ! In serial make sure the particle is in the domain
       ! unless we have the write_nan_outside override
       if (element<0 .and. .not.detector_list%write_nan_outside) then
          ewrite(-1,*) "Dealing with particle ", id, " named: ", trim(name)
          FLExit("Trying to initialise particle outside of computational domain")
       end if
    end if
    ! Otherwise, allocate and insert particle
    allocate(detector)
    allocate(detector%position(xfield%dim))
    allocate(detector%local_coords(local_coord_count(shape)))
    call insert(detector, detector_list)
    ! Populate particle
    detector%name=name
    detector%position=position
    detector%element=element
    detector%local_coords=lcoords
    detector%type=type
    detector%id_number=id

    allocate(detector%attributes(attribute_size(1)))
    allocate(detector%old_attributes(attribute_size(2)))
    allocate(detector%old_fields(attribute_size(3)))

    if (present(attribute_vals)) then
       detector%attributes = attribute_vals(1,1:attribute_size(1))
       detector%old_attributes = attribute_vals(2,1:attribute_size(2))
       detector%old_fields = attribute_vals(3,1:attribute_size(3))
    else
       detector%attributes(:) = 0
       detector%old_attributes(:) = 0
       detector%old_fields(:) = 0
    end if
    
  end subroutine create_single_particle

  subroutine move_particles(state, dt, timestep)
    !!Routine to loop over particle arrays and call move_lagrangian_detectors
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: dt
    integer, intent(in) :: timestep

    integer, dimension(3) :: attribute_size
    integer :: particle_groups, list_counter, dim
    integer, dimension(:), allocatable :: particle_arrays
    integer :: i, m, k
    integer :: nscalar, nvector, ntensor
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    call get_option("/geometry/dimension",dim)

    !Number of old_fields stored on particles
    nscalar = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')
    nvector = option_count('/material_phase/vector_field/prognostic/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/vector_field/prescribed/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/vector_field/diagnostic/particles/include_in_particles/store_old_field')
    ntensor = option_count('/material_phase/tensor_field/prognostic/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/tensor_field/prescribed/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/tensor_field/diagnostic/particles/include_in_particles/store_old_field')

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do
    ewrite(2,*), "In move_particles"
    list_counter = 1

    !Loop over all particle arrays
    do i = 1, particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          attribute_size(1)=0

          !Get attribute_size
          if (have_option(trim(subgroup_path) // "/attributes")) then
             attribute_size(1)=option_count(trim(subgroup_path) // "/attributes/scalar_attribute")
             attribute_size(1)=attribute_size(1)+dim*option_count(trim(subgroup_path) // "/attributes/vector_attribute")
             attribute_size(1)=attribute_size(1)+(dim**2)*option_count(trim(subgroup_path) // "/attributes/tensor_attribute")
          end if
          attribute_size(2)=0
          attribute_size(3)=0

          !Get attribute size of old_fields and old_attributes
          do m = 1,option_count(trim(subgroup_path) // "/attributes/scalar_attribute")
             if (have_option(trim(subgroup_path) // "/attributes/scalar_attribute["//int2str(m-1)//"]/python_fields")) then!!!check if this numbering works
                attribute_size(3) = nscalar + dim*nvector + ntensor*(dim**2)
                if (have_option(trim(subgroup_path) // "/attributes/scalar_attribute["//int2str(m-1)//"]/python_fields/store_old_attribute")) then
                   attribute_size(2)=attribute_size(2)+1
                end if
             end if
          end do
          do m = 1,option_count(trim(subgroup_path) // "/attributes/vector_attribute")
             if (have_option(trim(subgroup_path) // "/attributes/vector_attribute["//int2str(m-1)//"]/python_fields")) then!!!check if this numbering works
                attribute_size(3) = nscalar + dim*nvector + ntensor*(dim**2)
                if (have_option(trim(subgroup_path) // "/attributes/vector_attribute["//int2str(m-1)//"]/python_fields/store_old_attribute")) then
                   attribute_size(2)=attribute_size(2)+dim
                end if
             end if
          end do
          do m = 1,option_count(trim(subgroup_path) // "/attributes/tensor_attribute")
             if (have_option(trim(subgroup_path) // "/attributes/tensor_attribute["//int2str(m-1)//"]/python_fields")) then!!!check if this numbering works
                attribute_size(3) = nscalar + dim*nvector + ntensor*(dim**2)
                if (have_option(trim(subgroup_path) // "/attributes/tensor_attribute["//int2str(m-1)//"]/python_fields/store_old_attribute")) then
                   attribute_size(2)=attribute_size(2)+(dim**2)
                end if
             end if
          end do
          call move_lagrangian_detectors(state, particle_lists(list_counter), dt, timestep, attribute_size)
          list_counter = list_counter + 1
       end do
    end do
    deallocate(particle_arrays)

  end subroutine move_particles

  subroutine update_particle_subgroup_attributes_and_fields(state, time, dt, subgroup_path, p_list)
    !!Routine to set particle attributes 
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    real, intent(in) :: dt
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    type(detector_linked_list), intent(in) :: p_list
    
    character(len=PYTHON_FUNC_LEN) :: func
    character(len=FIELD_NAME_LEN), allocatable, dimension(:) :: field_name
    type(detector_type), pointer :: particle

    real, allocatable, dimension(:,:) :: positions
    real, allocatable, dimension(:,:) :: attribute_array
    real, allocatable, dimension(:,:) :: old_attributes
    character, allocatable, dimension(:,:) :: old_att_names
    character(len = OPTION_PATH_LEN) :: old_name

    real :: constant
    real, allocatable, dimension(:) :: vconstant
    real, allocatable, dimension(:,:) :: tconstant
    integer :: j, nparticles, l, m, n, dim, i
    integer :: nscalar, nvector, ntensor
    integer, allocatable, dimension(:) :: store_old_att
    real, allocatable, dimension(:,:) :: lcoords
    integer, allocatable, dimension(:) :: ele
    integer, dimension(3) :: n_oldatts !number of old scalar, vector and tensor attributes

    !Check if this processor contains particles
    nparticles = p_list%length

    if (nparticles==0) then
       return
    end if

    call get_option("/geometry/dimension",dim)

    !Set parameters to calculate attributes
    particle => p_list%first
    allocate(positions(size(particle%position),nparticles))
    allocate(attribute_array(size(particle%attributes),nparticles))
    allocate(lcoords(size(particle%local_coords),nparticles))
    allocate(ele(nparticles))
    allocate(old_attributes(size(particle%old_attributes),nparticles))

    n_oldatts(1) = option_count(trim(subgroup_path) //"/attributes/scalar_attribute/python_fields/store_old_attribute")
    n_oldatts(2) = option_count(trim(subgroup_path) //"/attributes/vector_attribute/python_fields/store_old_attribute")
    n_oldatts(3) = option_count(trim(subgroup_path) //"/attributes/tensor_attribute/python_fields/store_old_attribute")
    allocate(old_att_names(FIELD_NAME_LEN,sum(n_oldatts)))

    particle => p_list%first
    do j = 1,nparticles
       positions(:,j) = particle%position
       lcoords(:,j) = particle%local_coords
       ele(j) = particle%element
       old_attributes(:,j) = particle%old_attributes
       particle => particle%next
    end do

    !get number of scalar, vector and tensor attributes
    nscalar = option_count(trim(subgroup_path) // "/attributes/scalar_attribute")
    nvector = option_count(trim(subgroup_path) // "/attributes/vector_attribute")
    ntensor = option_count(trim(subgroup_path) // "/attributes/tensor_attribute")
    
    l=1
    particle => p_list%first
    !set old scalar attribute names
    do n = 1,nscalar
       if (have_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
          call get_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/name', old_name)
          old_att_names(1,l) = 'O'
          old_att_names(2,l) = 'l'
          old_att_names(3,l) = 'd'
          do j = 4,len_trim(old_name)+3
             old_att_names(j,l)=old_name(j-3:j-3)
          end do
          old_att_names(j,l) = C_NULL_CHAR
          l=l+1
       end if
    end do
    !set old vector attribute names
    do n = 1,nvector
       if (have_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
          call get_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/name', old_name)
          old_att_names(1,l) = 'O'
          old_att_names(2,l) = 'l'
          old_att_names(3,l) = 'd'
          do j = 4,len_trim(old_name)+3
             old_att_names(j,l)=old_name(j-3:j-3)
          end do
          old_att_names(j,l) = C_NULL_CHAR
          l=l+1
       end if
    end do
    !set old tensor attribute names
    do n = 1,ntensor
       if (have_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
          call get_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/name', old_name)
          old_att_names(1,l) = 'O'
          old_att_names(2,l) = 'l'
          old_att_names(3,l) = 'd'
          do j = 4,len_trim(old_name)+3
             old_att_names(j,l)=old_name(j-3:j-3)
          end do
          old_att_names(j,l) = C_NULL_CHAR
          l=l+1
       end if
    end do

    !Update scalar attributes, then vector attributes, then tensor attributes
    do n = 1,nscalar
       if (have_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/constant')) then
          call get_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/constant', constant)
          attribute_array(n,:) = constant
       else if (have_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/python')) then
          call get_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/python', func)
          call set_particle_scalar_attribute_from_python(attribute_array(n,:), positions(:,:), nparticles, func, time, dt)
       else if (have_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/python_fields')) then
          call get_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/python_fields', func)
          call set_particle_scalar_attribute_from_python_fields(p_list, state, positions(:,:), lcoords(:,:), ele(:), nparticles, &
               & attribute_array(n,:), n_oldatts, old_att_names, old_attributes, func, time, dt)
       else if (have_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/from_checkpoint_file')) then
          particle => p_list%first
          attribute_array(n,:) = particle%attributes(n)
       end if
    end do
    allocate(vconstant(dim))
    do n = 1,nvector
       if (have_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/constant')) then
          call get_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/constant', vconstant)
          do i = 1, dim
             attribute_array(nscalar+i+((n-1)*dim),:) = vconstant(i)
          end do
       else if (have_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/python')) then
          call get_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/python', func)
          call set_particle_vector_attribute_from_python(attribute_array(nscalar+1+((n-1)*dim):nscalar+dim+((n-1)*dim),:), positions(:,:), nparticles, func, time, dt)
       else if (have_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/python_fields')) then
          call get_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/python_fields', func)
          call set_particle_vector_attribute_from_python_fields(p_list, state, positions(:,:), lcoords(:,:), ele(:), nparticles, &
               & attribute_array(nscalar+1+((n-1)*dim):nscalar+dim+((n-1)*dim),:), n_oldatts, old_att_names, old_attributes, func, time, dt)
       else if (have_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/from_checkpoint_file')) then
          particle => p_list%first
          do i = 1, dim
             attribute_array(nscalar+i+((n-1)*dim),:) = particle%attributes(nscalar+i+((n-1)*dim))
          end do
       end if
    end do
    allocate(tconstant(dim,dim))
    do n = 1,ntensor
       if (have_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/constant')) then
          call get_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/constant', tconstant)
          do i = 1,dim
             do j = 1,dim
                attribute_array(nscalar+(nvector*dim)+((i-1)*dim)+j+((n-1)*(dim**2)),:) = tconstant(j,i)
             end do
          end do
       else if (have_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/python')) then
          call get_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/python', func)
          call set_particle_tensor_attribute_from_python(attribute_array(nscalar+(nvector*dim)+1+((n-1)*(dim**2)):nscalar+(nvector*dim)+(dim**2)+((n-1)*(dim**2)),:), &
               positions(:,:), nparticles, func, time, dt)
       else if (have_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/python_fields')) then
          call get_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/python_fields', func)
          call set_particle_tensor_attribute_from_python_fields(p_list, state, positions(:,:), lcoords(:,:), ele(:), nparticles, &
               & attribute_array(nscalar+(nvector*dim)+1+((n-1)*(dim**2)):nscalar+(nvector*dim)+(dim**2)+((n-1)*(dim**2)),:), &
               n_oldatts, old_att_names, old_attributes, func, time, dt)
       else if (have_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/from_checkpoint_file')) then
          particle => p_list%first
          do i = 1,(dim**2)
             attribute_array(nscalar+(nvector*dim)+i+((n-1)*(dim**2)),:) = particle%attributes(nscalar+(nvector*dim)+i+((n-1)*(dim**2)))
          end do
       end if
    end do
    
    !Set attribute values and old_attribute values
    particle => p_list%first
    if (size(particle%old_attributes)==0) then
       do j = 1,nparticles
          particle%attributes = attribute_array(:,j)
          particle => particle%next
       end do
    else
       allocate(store_old_att(size(particle%attributes)))
       !copy old scalar attributes
       do n = 1,nscalar
          if (have_option(trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
             store_old_att(n)=1
          else
             store_old_att(n)=0
          end if
       end do
       !copy old vector attributes
       do n = 1,nvector
          if (have_option(trim(subgroup_path) // '/attributes/vector_attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
             store_old_att(nscalar+1+((n-1)*dim):nscalar+dim+((n-1)*dim))=1
          else
             store_old_att(nscalar+1+((n-1)*dim):nscalar+dim+((n-1)*dim))=0
          end if
       end do
       !copy old tensor attributes
       do n = 1,ntensor
          if (have_option(trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
             store_old_att(nscalar+(nvector*dim)+1+((n-1)*(dim**2)):nscalar+(nvector*dim)+(dim**2)+((n-1)*(dim**2)))=1
          else
             store_old_att(nscalar+(nvector*dim)+1+((n-1)*(dim**2)):nscalar+(nvector*dim)+(dim**2)+((n-1)*(dim**2)))=0
          end if
       end do
       
       do j = 1,nparticles
          particle%attributes = attribute_array(:,j)
          m=1
          do n = 1,size(particle%attributes)
             if (store_old_att(n)==1) then
                particle%old_attributes(m) = particle%attributes(n)
                m=m+1
             end if
          end do
          particle => particle%next
       end do
    end if

    particle => p_list%first
    if (size(particle%old_fields)/=0) then
       call update_particle_subgroup_fields(state, ele, lcoords, p_list)
    end if

    deallocate(positions)
    deallocate(lcoords)
    deallocate(ele)
    deallocate(attribute_array)
    deallocate(old_attributes)
    deallocate(old_att_names)

  end subroutine update_particle_subgroup_attributes_and_fields

  subroutine update_particle_attributes_and_fields(state, time)
    !!Routine to loop over particle arrays and update particle attributes
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    real :: dt
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: particle
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    integer :: i, k
    integer :: particle_groups, list_counter
    integer, dimension(:), allocatable :: particle_arrays

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    ewrite(2,*), "In update_particle_attributes_and_fields"

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/timestepping/timestep", dt)
    list_counter = 1

    !Update particle attributes by array
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          if (particle_lists(list_counter)%length==0) then
             list_counter = list_counter + 1
             cycle
          end if
          particle => particle_lists(list_counter)%first
          if (size(particle%attributes)/=0) then
             call update_particle_subgroup_attributes_and_fields(state, time, dt, subgroup_path, particle_lists(list_counter))
          end if
          list_counter = list_counter + 1
       end do
    end do
  end subroutine update_particle_attributes_and_fields

  subroutine update_particle_subgroup_fields(state, ele, lcoords, p_list)
    
    type(state_type), dimension(:), intent(in) :: state
    real, dimension(:,:), intent(in) :: lcoords
    integer, dimension(:), intent(in) :: ele
    type(detector_linked_list), intent(in) :: p_list

    character(len = OPTION_PATH_LEN) :: name
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    real, allocatable, dimension(:,:) :: old_field_vals
    real :: value
    real, allocatable, dimension(:) :: vvalue
    real, allocatable, dimension(:,:) :: tvalue
    type(detector_type), pointer :: particle
    integer :: phase, f, l, j, dim, i, k
    integer :: s_oldfield, v_oldfield, t_oldfield

    call get_option("/geometry/dimension",dim)
    
    particle => p_list%first
    allocate(old_field_vals(size(particle%old_fields),size(lcoords(1,:))))!particle old fields and number of particles
    allocate(vvalue(dim))
    allocate(tvalue(dim,dim))
    do phase=1,size(state)
       l=1
       !Scalar old fields
       do f = 1, size(state(phase)%scalar_names)
          sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(f))
          if (sfield%option_path=="".or.aliased(sfield)) then
             cycle
          else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
             do j = 1,size(lcoords(1,:))
                value = eval_field(ele(j), sfield, lcoords(:,j))
                old_field_vals(l,j)=value
             end do
             l=l+1
          end if
       end do
    end do

    do phase=1,size(state)
       s_oldfield = l-1
       l=1
       !Vector old fields
       do f = 1, size(state(phase)%vector_names)
          vfield => extract_vector_field(state(phase),state(phase)%vector_names(f))
          if (vfield%option_path=="".or.aliased(vfield)) then
             cycle
          else if (have_option(trim(complete_field_path(vfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
             do j = 1,size(lcoords(1,:))
                vvalue(:) = eval_field(ele(j), vfield, lcoords(:,j))
                old_field_vals(s_oldfield+1+((l-1)*dim):s_oldfield+dim+((l-1)*dim),j)=vvalue(:)
             end do
             l=l+1
          end if
       end do
    end do

    do phase=1,size(state)
       v_oldfield = l-1
       l=1
       !Tensor old fields
       if (associated(state(phase)%tensor_names)) then
          do f = 1, size(state(phase)%tensor_names)
             tfield => extract_tensor_field(state(phase),state(phase)%tensor_names(f))
             if (tfield%option_path=="".or.aliased(tfield)) then
                cycle
             else if (have_option(trim(complete_field_path(tfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
                do j = 1,size(lcoords(1,:))
                   tvalue(:,:) = eval_field(ele(j), tfield, lcoords(:,j))
                   do k = 1,dim
                      old_field_vals(s_oldfield+(v_oldfield*dim)+((l-1)*(dim**2))+((k-1)*dim)+1:s_oldfield+(v_oldfield*dim)+((l-1)*(dim**2))+((k-1)*dim)+dim,j)=tvalue(:,k)!!!check this
                   end do
                end do
                l=l+1
             end if
          end do
       end if
    end do

    particle => p_list%first
    do j = 1,size(lcoords(1,:))
       particle%old_fields=old_field_vals(:,j)
       particle=>particle%next
    end do

    deallocate(old_field_vals)
    deallocate(vvalue)
    deallocate(tvalue)

  end subroutine update_particle_subgroup_fields

  subroutine write_particles_loop(state, time, dt)
    !!Subroutine to loop over particle_lists and call write_particles for each list
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt

    integer :: attribute_dims
    integer :: i, k, dim
    integer :: particle_groups, list_counter
    integer :: s_att, v_att, t_att, tot_atts
    integer, dimension(:), allocatable :: particle_arrays
    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    call get_option("/geometry/dimension",dim)

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    ewrite(1,*) "In write_particles_loop"
    
    list_counter = 1
    do i = 1, particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          s_att = option_count(trim(subgroup_path) // '/attributes/scalar_attribute')
          v_att = option_count(trim(subgroup_path) // '/attributes/vector_attribute')
          t_att = option_count(trim(subgroup_path) // '/attributes/tensor_attribute')
          attribute_dims=s_att + dim*v_att + (dim**2)*t_att
          call write_particles_subgroup(state, particle_lists(list_counter), attribute_dims, time, dt)
          list_counter = list_counter + 1
       end do
    end do

    deallocate(particle_arrays)

  end subroutine write_particles_loop

  subroutine write_particles_subgroup(state, detector_list, attribute_dims, time, dt)
    !!< Write values of particles to the previously opened particles file.
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: time, dt
    integer, intent(in) :: attribute_dims !dimensions of particles attribute information carried (attributes at current timestep, field values and attribute values at previous timestep)

    character(len=10) :: format_buffer
    integer :: i, check_no_det, totaldet_global
    type(detector_type), pointer :: detector

    ewrite(1,*) "In write_particles"

    !Computing the global number of particles. This is to prevent hanging
    !when there are no particles on any processor
    check_no_det=1
    if (detector_list%length==0) then
       check_no_det=0
    end if
    call allmax(check_no_det)
    if (check_no_det==0) then
       return
    end if

    ! If binary output use this:
    if (detector_list%binary_output) then    
       call write_mpi_out_particles(state,detector_list, attribute_dims, time,dt)
    else ! This is only for single processor with ascii output
       if(getprocno() == 1) then
          format_buffer=reals_format(1)
          write(detector_list%output_unit, format_buffer, advance="no") time
          write(detector_list%output_unit, format_buffer, advance="no") dt
       end if

       ! Next columns contain the positions of all the particles.
       detector => detector_list%first
       positionloop: do i=1, detector_list%length
          format_buffer=reals_format(size(detector%position))
          write(detector_list%output_unit, format_buffer, advance="no") &
               detector%position
          detector => detector%next
       end do positionloop

       ! Next columns contain the attributes of particles
       detector => detector_list%first
       attributeloop: do i=1,detector_list%length
          if (attribute_dims/=0) then
             format_buffer=reals_format(attribute_dims)
             write(detector_list%output_unit, format_buffer, advance="no") &
                  detector%attributes
          end if
          detector => detector%next
       end do attributeloop
       
       ! Output end of line
       flush(detector_list%output_unit)
    end if

    totaldet_global=detector_list%length
    call allsum(totaldet_global)
    ewrite(2,*) "Found", detector_list%length, "local and", totaldet_global, "global detectors"
    
    if (totaldet_global/=detector_list%total_num_det) then
       ewrite(2,*) "We have either duplication or have lost some det"
       ewrite(2,*) "totaldet_global", totaldet_global
       ewrite(2,*) "total_num_det", detector_list%total_num_det
    end if
    
    ewrite(1,*) "Exiting write_particles"
    
  contains
    
    function reals_format(reals)
      character(len=10) :: reals_format
      integer :: reals
      
      write(reals_format, '(a,i0,a)') '(',reals,'e15.6e3)'
      
    end function reals_format
    
  end subroutine write_particles_subgroup

  subroutine write_mpi_out_particles(state, detector_list, attribute_dims, time, dt)
    !!< Writes particle information (position, attributes, etc.) into particle file using MPI output 
    ! commands so that when running in parallel all processors can write at the same time information into the file at the right location.       

    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: time, dt
    integer, intent(in) :: attribute_dims

    integer :: i, ierror, realsize, dim, procno
    integer(KIND = MPI_OFFSET_KIND) :: location_to_write, offset
    integer :: number_total_columns

    type(vector_field), pointer :: vfield
    type(detector_type), pointer :: node

    ewrite(2, *) "In write_mpi_out_particles"

    detector_list%mpi_write_count = detector_list%mpi_write_count + 1
    ewrite(2, *) "Writing particle output ", detector_list%mpi_write_count
    
    procno = getprocno()


    call mpi_type_extent(getpreal(), realsize, ierror)
    assert(ierror == MPI_SUCCESS)

    vfield => extract_vector_field(state, "Coordinate")
    dim = vfield%dim
                           ! Time data
    number_total_columns = 2 + &
                           ! Particle coordinates
                         & detector_list%total_num_det * dim + &
                           ! Particle attributes
                         & detector_list%total_num_det * attribute_dims

    ! raise kind of one of the variables (each individually is a 4 byte-integer) such that the calculation is coerced to be of MPI_OFFSET_KIND (typically 8 bytes)
    ! this is necessary for files bigger than 2GB
    location_to_write = (int(detector_list%mpi_write_count, kind=MPI_OFFSET_KIND) - 1) * number_total_columns * realsize

    if(procno == 1) then
      ! Output time data
      call mpi_file_write_at(detector_list%mpi_fh, location_to_write, time, 1, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
        
      call mpi_file_write_at(detector_list%mpi_fh, location_to_write + realsize, dt, 1, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
    end if
    location_to_write = location_to_write + 2 * realsize

    node => detector_list%first
    position_loop: do i = 1, detector_list%length
      ! Output detector coordinates
      assert(size(node%position) == dim)  
    
      offset = location_to_write + (node%id_number - 1) * dim * realsize

      call mpi_file_write_at(detector_list%mpi_fh, offset, node%position, dim, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
      node => node%next
    end do position_loop
    assert(.not. associated(node))
    location_to_write = location_to_write + detector_list%total_num_det * dim * realsize

    if (attribute_dims/=0) then
       node => detector_list%first
       attribute_loop: do i = 1, detector_list%length
          assert(size(node%attributes) == attribute_dims)
          
          offset = location_to_write + (node%id_number - 1) * attribute_dims * realsize
          call mpi_file_write_at(detector_list%mpi_fh, offset, node%attributes, attribute_dims, getpreal(), MPI_STATUS_IGNORE, ierror)
          assert(ierror == MPI_SUCCESS)
          node => node%next
       end do attribute_loop
       assert(.not. associated(node))
    end if
    
    call mpi_file_sync(detector_list%mpi_fh, ierror)
    assert(ierror == MPI_SUCCESS)

    ewrite(2, *) "Exiting write_mpi_out"
   
  end subroutine write_mpi_out_particles

  subroutine checkpoint_particles_loop(state,prefix,postfix,cp_no)
    !!Subroutine to loop over particle_lists and call checkpoint_particles for each list
    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), intent(in) :: postfix
    integer, optional, intent(in) :: cp_no !Checkpoint number of the simulation
    
    character(len = OPTION_PATH_LEN) :: lpostfix
    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path, subgroup_path_name, name

    integer, dimension(3) :: attribute_size
    integer :: nscalar, nvector, ntensor
    integer :: s_att, v_att, t_att, tot_atts
    integer, dimension(:), allocatable :: particle_arrays
    integer :: i, m, k, particle_groups, list_counter, dim

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    call get_option("/geometry/dimension",dim)

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    ewrite(1, *) "Checkpointing particles"

    assert(len_trim(prefix) > 0)

    lpostfix = postfix

    !Number of old_fields stored on particles
    nscalar = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')
    nvector = option_count('/material_phase/vector_field/prognostic/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/vector_field/prescribed/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/vector_field/diagnostic/particles/include_in_particles/store_old_field')
    ntensor = option_count('/material_phase/tensor_field/prognostic/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/tensor_field/prescribed/particles/include_in_particles/store_old_field') &
    + option_count('/material_phase/tensor_field/diagnostic/particles/include_in_particles/store_old_field')

    list_counter = 1
    !Loop over all particle arrays
    do i = 1, particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          subgroup_path_name = trim(group_path) // "/particle_subgroup::"!!option path used in update_particle_options
          attribute_size(1)=0
          !Get attribute_size
          s_att = option_count(trim(subgroup_path) // '/attributes/scalar_attribute')
          v_att = option_count(trim(subgroup_path) // '/attributes/vector_attribute')
          t_att = option_count(trim(subgroup_path) // '/attributes/tensor_attribute')
          if (s_att + v_att + t_att > 0) then
             attribute_size(1)= s_att + dim*v_att + (dim**2)*t_att
          end if
          attribute_size(2)=0
          attribute_size(3)=0
          !Get attribute size of old_fields and old_attributes
          do m = 1,s_att
             if (have_option(trim(subgroup_path) // "/attributes/scalar_attribute["//int2str(m-1)//"]/python_fields")) then!!!check if this numbering works
                attribute_size(3) = nscalar + dim*nvector + ntensor*(dim**2)
                if (have_option(trim(subgroup_path) // "/attributes/scalar_attribute["//int2str(m-1)//"]/python_fields/store_old_attribute")) then
                   attribute_size(2)=attribute_size(2)+1
                end if
             end if
          end do
          do m = 1,v_att
             if (have_option(trim(subgroup_path) // "/attributes/vector_attribute["//int2str(m-1)//"]/python_fields")) then!!!check if this numbering works
                attribute_size(3) = nscalar + dim*nvector + ntensor*(dim**2)
                if (have_option(trim(subgroup_path) // "/attributes/vector_attribute["//int2str(m-1)//"]/python_fields/store_old_attribute")) then
                   attribute_size(2)=attribute_size(2)+dim
                end if
             end if
          end do
          do m = 1,t_att
             if (have_option(trim(subgroup_path) // "/attributes/tensor_attribute["//int2str(m-1)//"]/python_fields")) then!!!check if this numbering works
                attribute_size(3) = nscalar + dim*nvector + ntensor*(dim**2)
                if (have_option(trim(subgroup_path) // "/attributes/tensor_attribute["//int2str(m-1)//"]/python_fields/store_old_attribute")) then
                   attribute_size(2)=attribute_size(2)+(dim**2)
                end if
             end if
          end do
          call get_option(trim(subgroup_path) // "/name", name)
          tot_atts = s_att + v_att + t_att
          call checkpoint_particles_subgroup(state,prefix,lpostfix,cp_no,particle_lists(list_counter),attribute_size, tot_atts, name, subgroup_path_name)
          list_counter = list_counter + 1
       end do
    end do

    deallocate(particle_arrays)

  end subroutine checkpoint_particles_loop

  subroutine checkpoint_particles_subgroup(state,prefix,lpostfix,cp_no,particle_list,attribute_size, tot_atts, name, subgroup_path_name)
    !!<Checkpoint Particles

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), intent(in) :: lpostfix
    integer, optional, intent(in) :: cp_no !Checkpoint number of the simulation
    character(len = *), intent(in) :: name
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path_name

    type(detector_linked_list), intent(inout) :: particle_list
    integer, dimension(3), intent(in) :: attribute_size
    integer, intent(in) :: tot_atts

    integer(KIND=MPI_OFFSET_KIND) :: location_to_write, offset
    
    type(detector_type), pointer :: node
    character(len = OPTION_PATH_LEN) :: particles_cp_filename
    type(vector_field), pointer :: vfield

    integer, ALLOCATABLE, DIMENSION(:) :: status
    real, dimension(:), allocatable :: buffer
    
    integer, save :: fhdet=0
    integer :: j, ierror
    integer :: nints, realsize, dimen, num_particles, number_total_columns, intsize

    integer :: num_procs, proc_num

    num_procs = getnprocs()
    proc_num = getprocno()

    !num_particles = particle_list%total_num_det
    num_particles = particle_list%length

    ! Construct a new particle checkpoint filename
    !get name of particle array here to construct the output file
    particles_cp_filename = trim(prefix)
    if(present(cp_no)) particles_cp_filename = trim(particles_cp_filename) // "_" // int2str(cp_no)
    particles_cp_filename = trim(particles_cp_filename) // "_" // trim(lpostfix)

    !!< Writes particle last position into particles file using MPI output 
    ! commands so that when running in parallel all processors can write at the same time information into the file at the right location.
    
    call MPI_FILE_OPEN(MPI_COMM_SELF, trim(particles_cp_filename) // '_particles.' // trim(parallel_filename(name)) // '.attributes.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, fhdet, ierror)

    ewrite(1,*) "after opening the ierror is:", ierror

    allocate( status(MPI_STATUS_SIZE) )

    call MPI_TYPE_EXTENT(getpinteger(), intsize, ierror)
    call MPI_TYPE_EXTENT(getpreal(), realsize, ierror)

    vfield => extract_vector_field(state(1),"Velocity")

    dimen=vfield%dim

    number_total_columns=num_particles*(dimen+sum(attribute_size))

    node => particle_list%first

    offset=0*intsize

    call MPI_FILE_WRITE_AT(fhdet,offset,particle_list%length,1,getpinteger(),status,IERROR)

    location_to_write=1*intsize

    positionloop_cp: do j=1, particle_list%length
      offset = location_to_write+(j-1)*(size(node%position)+sum(attribute_size))*realsize
      ewrite(1,*) "after file set view position ierror is:", ierror

      allocate(buffer(size(node%position)+sum(attribute_size)))
      buffer(1:size(node%position))=node%position
      if (attribute_size(1)/=0) then
         buffer(1+size(node%position):size(node%position)+attribute_size(1))=node%attributes
      end if
      if (attribute_size(2)/=0) then
         buffer(1+size(node%position)+attribute_size(1):size(node%position)+attribute_size(1) &
              +attribute_size(2))=node%old_attributes
      end if
      if (attribute_size(3)/=0) then
         buffer(1+size(node%position)+attribute_size(1)+attribute_size(2):size(node%position)+attribute_size(1) &
              +attribute_size(2)+attribute_size(3))=node%old_fields
      end if
      nints=size(node%position)+sum(attribute_size)
      
      call MPI_FILE_WRITE_AT(fhdet,offset,buffer,nints,getpreal(),status,ierror)

      ewrite(1,*) "after sync position ierror is:", ierror
      deallocate(buffer)
      node => node%next
    end do positionloop_cp

    call update_particle_subgroup_options(trim(particles_cp_filename) // "_particles", "binary", particle_list, name, tot_atts, subgroup_path_name)
    
    if (fhdet/=0) then
       call MPI_FILE_CLOSE(fhdet, ierror)
       if (ierror/=0) then
          ewrite(0,*) "Warning: failed to close .particles checkpoint file open with mpi_file_open"
       end if
    end if
    
  end subroutine checkpoint_particles_subgroup

  subroutine update_particle_subgroup_options(filename, format, particle_list, name, tot_atts, subgroup_path_name)
    !! Updates the initial options of particles in the schema file for reinitialization after checkpointing.
    !! Updates schema options for the initial number of particles and their initial positions. 

    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: format
    character(len = *), intent(in) :: name

    type(detector_linked_list), intent(inout) :: particle_list
    integer, intent(in) :: tot_atts
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path_name

    integer :: num_particles, j, stat
    logical :: particles_s, particles_v, particles_t

    character(len = 254) :: temp_string

    num_particles = particle_list%total_num_det

    temp_string=name

    ewrite(1,*) 'In update_particles_options'
    ewrite(1,*) temp_string
    
    call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/number_of_particles")
    call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/initial_position")
    
    call set_option(trim(subgroup_path_name) // trim(temp_string) // "/number_of_particles/", &
         & num_particles, stat = stat)
    
    assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))
    
    call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/initial_position/from_file/file_name", trim(filename)// "." // trim(temp_string), stat)
    
    if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
       FLAbort("Failed to set particles options filename when checkpointing particles with option path " // "/particles/particle_array::" // trim(temp_string))
    end if
    
    call set_option(trim(subgroup_path_name) // trim(temp_string) // "/initial_position/from_file/format/", trim(format), stat)
    
    if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
       FLAbort("Failed to set particles options format when checkpointing particles with option path " // "/particles/particle_group")
    end if
    
    do j = 1, tot_atts
       particles_s = have_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/scalar_attribute["//int2str(j-1)//"]/constant")
       particles_v = have_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/vector_attribute["//int2str(j-1)//"]/constant")
       particles_t = have_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/tensor_attribute["//int2str(j-1)//"]/constant")
       if (particles_s) then
          call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/scalar_attribute["//int2str(j-1)//"]/constant")
          call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/attributes/scalar_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
             FLAbort("Failed to set scalar field particles filename when checkpointing")
          end if
          call set_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/scalar_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/format/", trim(format), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
             FLAbort("Failed to set scalar field particles options format when checkpointing")
          end if
       else if (particles_v) then
          call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/vector_attribute["//int2str(j-1)//"]/constant")
          call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/attributes/vector_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
             FLAbort("Failed to set vector field particles filename when checkpointing")
          end if
          call set_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/vector_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/format/", trim(format), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
             FLAbort("Failed to set vector field particles options format when checkpointing")
          end if
       else if (particles_t) then
          call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/tensor_attribute["//int2str(j-1)//"]/constant")
          call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/attributes/tensor_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
             FLAbort("Failed to set tensor field particles filename when checkpointing")
          end if
          call set_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/tensor_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/format/", trim(format), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
             FLAbort("Failed to set tensor field particles options format when checkpointing")
          end if
       end if
    end do
    
  end subroutine update_particle_subgroup_options

  subroutine get_particles(p_array, p_allocated)
    !Send particle arrays to another routine

    type(detector_linked_list), allocatable, dimension(:), intent(out) :: p_array
    integer, intent(out) :: p_allocated

    integer :: i, particle_groups

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) then
       FLAbort("No particle groups exist")
       return
    end if
    
    if (allocated(particle_lists)) then
       p_allocated = 1
       allocate(p_array(size(particle_lists)))
       do i = 1,size(particle_lists)
          p_array(i) = particle_lists(i)
       end do
    else
       p_allocated = 0
    end if
    
  end subroutine get_particles

  subroutine get_particle_arrays(lgroup, group_arrays, group_attribute, lattribute)
    !Read in a particle group and attribute name or particle subgroup, send back numbers of particle arrays and particle attribute

    character(len=OPTION_PATH_LEN), intent(in) :: lgroup
    character(len=OPTION_PATH_LEN), optional, intent(in) :: lattribute
    integer, allocatable, dimension(:), intent(out) :: group_arrays
    integer, optional, intent(out) :: group_attribute

    character(len=OPTION_PATH_LEN) :: group_name, attribute_name, subgroup_name
    integer :: particle_groups, array_counter, particle_subgroups, particle_attributes
    integer :: i, j, k, l
    

    logical :: found_attribute
    
    particle_groups = option_count("/particles/particle_group")

    found_attribute = .false.
    array_counter = 0
    do i = 1, particle_groups
       call get_option("/particles/particle_group["//int2str(i-1)//"]/name", group_name)
       particle_subgroups = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
       if (trim(group_name)==trim(lgroup)) then
          allocate(group_arrays(particle_subgroups))
          if (present(lattribute)) then
             particle_attributes = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)//"]/attributes/scalar_attribute")
             do k = 1, particle_attributes
                call get_option("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)// &
                     "]/attributes/scalar_attribute["//int2str(k-1)//"]/name", attribute_name)
                if (trim(attribute_name)==trim(lattribute)) then
                   found_attribute = .true.
                   group_attribute = k
                end if
             end do
             if (found_attribute.eqv..false.) then
                FLExit("Could not find particle attribute "//trim(lattribute)//" in particle group "//trim(lgroup)//". Check attribute is a scalar.")
             end if
          end if
          j=1
          do l = array_counter+1, array_counter+particle_subgroups
             group_arrays(j) = l
             j=j+1
          end do
          return
       end if
       array_counter = array_counter + particle_subgroups
    end do
    
    FLExit("Could not find particle group "//trim(lgroup))
    
    
  end subroutine get_particle_arrays

  subroutine update_list_lengths(list_num)
    integer, intent(in) :: list_num

    particle_lists(list_num)%total_num_det = particle_lists(list_num)%total_num_det + 1
    particle_lists(list_num)%length = particle_lists(list_num)%length + 1

  end subroutine update_list_lengths
  subroutine destroy_particles()
    type(detector_linked_list), pointer :: del_particle_lists
    integer :: i
    !Deallocate all particle arrays (detector lists)
    if (allocated(particle_lists)) then
       do i=1,size(particle_lists)
          del_particle_lists => particle_lists(i)
          call deallocate(del_particle_lists)
       end do
       deallocate(particle_lists)
    end if
    
  end subroutine destroy_particles

end module particles

  
    
    
    
