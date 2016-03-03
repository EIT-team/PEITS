#ifndef ELECTRODE_HELPERS_HH
#define ELECTRODE_HELPERS_HH

using namespace std;

struct CurrentProtocol
{
  // constructor: loads current protocol from file
  CurrentProtocol(const char *filename)
  {
    error_ = 0;
    FILE *F = 0;
    F = fopen(filename, "r");
    if (F == NULL)
    {
      cout << "Protocol file could not be loaded!" << endl;
      error_ = 1;
      return;
    }

	  while(!std::feof(F))
	  {
		  int d1,d2,m1,m2;
		  int error = fscanf(F, "%d,%d,%d,%d\n", &d1, &d2, &m1, &m2);
		  if (!error)
			  cout << "File could not be read: " << F << endl;
		  vector<int> measurement(4);
		  measurement[0]=d1;
		  measurement[1]=d2;
		  measurement[2]=m1;
		  measurement[3]=m2;
		  current_protocol_.push_back(measurement);
	  }
	  fclose(F);

    uniquePatterns();
	  input_current_ = Dune::Fem::Parameter::getValue< double>( "input.current" );
  }

  void uniquePatterns()
  /* this function calculates all unique drive and measurement patterns
     and stores them in unique_patterns_. In retrace_unique_patterns_
     we can then find the corresponding solution for a step in the current
     protocol and scale it with current (1 if measurement and current if drive) */
  {
	  vector<int> pair(2), invpair(2), retrace(2);
    pair[0] = current_protocol_[0][0];
	  pair[1] = current_protocol_[0][1];
	  unique_patterns_.push_back(pair); retrace[0] = 0;
    pair[0] = current_protocol_[0][2];
	  pair[1] = current_protocol_[0][3];
    reverse_copy(pair.begin(),pair.end(),invpair.begin());
    if( pair == unique_patterns_[0] )
		  {retrace[1] = 0;}
    else if (invpair == unique_patterns_[0])
      {retrace[1] = -1;}
    else
	    {unique_patterns_.push_back(pair); retrace[1] = 1;}
    retrace_current_pattern_.push_back(retrace);  // this can be optimised: if it is the same than drive, do not write it

	  for(unsigned int i = 1 ; i<current_protocol_.size() ; ++i)
	  {
	    pair[0] = current_protocol_[i][0];
	    pair[1] = current_protocol_[i][1];
      reverse_copy(pair.begin(),pair.end(),invpair.begin());
	    for(unsigned int j = 0 ; j<unique_patterns_.size() ; ++j)
	    {
	      if( pair == unique_patterns_[j] )
		      {retrace[0] = j; break;}
        else if (invpair == unique_patterns_[j])
          {retrace[0] = -j-1; break;}  // code opposite injection by a minus sign. for pattern 0 to be unique, we have to shift the negative ones by 1
	      else if( j == unique_patterns_.size() -1 )
		      {unique_patterns_.push_back(pair); retrace[0] = j + 1;}
	    }

      pair[0] = current_protocol_[i][2];
	    pair[1] = current_protocol_[i][3];
      reverse_copy(pair.begin(),pair.end(),invpair.begin());
	    for(unsigned int j = 0 ; j<unique_patterns_.size() ; ++j)
	    {
	      if( pair == unique_patterns_[j] )
		      {retrace[1] = j; break;}
        else if (invpair == unique_patterns_[j])
          {retrace[1] = -j-1; break;}
	      else if( j == unique_patterns_.size() -1 )
		      {unique_patterns_.push_back(pair); retrace[1] = j + 1;}
	    }
      retrace_current_pattern_.push_back(retrace);
	  }
  }

  
  double current( int electrodenumber, int pattern_no )
  {
	  if (electrodenumber == unique_patterns_[pattern_no][0])
		  return 1.0;
		else if (electrodenumber == unique_patterns_[pattern_no][1])
		  return -1.0;
		else
		  return 0;
  }

  vector< vector<int> > current_protocol_;
  vector< vector<int> > unique_patterns_;
  vector< vector<int> > retrace_current_pattern_;
  double input_current_;
  int error_;
};


struct ElectrodePositions
{
  // constructor: loads electrode positions from file
  ElectrodePositions(const char *filename)
  {
    error_ = 0;
	  FILE *F = 0;
	  F = fopen(filename, "r");
    if (F == NULL)
    {
      cout << "Electrode positions file could not be loaded!" << endl;
      error_ = 1;
      return;
    }

	  while(!feof(F))
	  {

      // when electrodes are defined by nodes, these nodes are stored
      if ( Dune::Fem::Parameter::getValue< bool >("electrode.use_node_assignment") )
      {
        vector<double> bla;
        electrode_positions_.push_back(bla);
        double node;
        char delimiter;
        // read line by line and store as double => for later use these node indices have to be converted to int
        while ( fscanf(F, "%lf%c", &node, &delimiter) != EOF )
        {
          electrode_positions_.back().push_back(node-1); // matlab starts counting at 1 and dune at 0
          if (delimiter == '\n')
            break;
        }
        if (electrode_positions_.back().size() == 0) // this is necessary because Matlab writes a \n after last entry.
        {
          electrode_positions_.pop_back();
          break;
        }
        const double contactimpedance = Dune::Fem::Parameter::getValue< double >( "contact.impedance" );
        contactimpedance_.push_back(contactimpedance);
      }
      else // otherwise the electrode positions are written
      {
		    double x,y,z;
		    int error = fscanf(F, "%lf,%lf,%lf\n", &x, &y, &z);
		    if (!error)
		      cout << "File could not be read: " << filename << endl;
		    vector<double> center(3);
		    center[0]=x;
		    center[1]=y;
		    center[2]=z;
		    electrode_positions_.push_back(center);
        const double contactimpedance = Dune::Fem::Parameter::getValue< double >( "contact.impedance" );
        contactimpedance_.push_back(contactimpedance);
      }

	  }
	  fclose(F);

    const string diameterfilename = Dune::Fem::Parameter::getValue< string >( "electrode.diameter_file" );
    if (diameterfilename.compare("none") != 0)
    {
      cout << "using individual diameter specifications for electrodes" << endl;
      F = fopen(diameterfilename.c_str(), "r");
      if (F == NULL)
      {
        cout << "Electrode diameter file could not be loaded!" << endl;
        error_ = 1;
        return;
      }
      while(!feof(F))
      {
        double x;
        int error = fscanf(F, "%lf\n", &x);
        if (!error)
          cout << "File could not be read: " << diameterfilename << endl;
        electrode_diameter_.push_back(x);
      }
      fclose(F);
      if (electrode_diameter_.size() != electrode_positions_.size())
        cout << "The number of electrodes is different in diameter file!!" << endl;
    }
    else
    {
      cout << "using the same diameter for all electrodes" << endl;
      for ( unsigned int bla = 0; bla < electrode_positions_.size(); bla++)
        electrode_diameter_.push_back(Dune::Fem::Parameter::getValue< double >( "electrode.diameter" ));
    }
  }

  template< class FieldVector >
  int electrodeNumber (const FieldVector &center) const
  {
	  int number = 0;
	  for ( unsigned int i=0 ; i<electrode_positions_.size() ; ++i)
	  {
		  if ( onElectrode(center,i+1) )
		  {
			  number = i+1;
			  continue;
		  }
	  }
	  return number;
  }


  // if electrodes are defined by nodes, then use this version
  int electrodeNumber_node (const vector<int> &corners) const
  {
    int number = 0;
	  for ( unsigned int i=0 ; i<electrode_positions_.size() ; ++i)
	  {
		  if ( onElectrode_node(corners,i+1) )
		  {
			  number = i+1;
			  continue;
		  }
	  }
	  return number;
  }


  template< class FieldVector >
  int fixToElectrodeNumber (const FieldVector &center)
  {
	  int number = 0;
	  for ( unsigned int i=0 ; i<electrode_positions_.size() ; ++i)
	  {
		  if ( fixToElectrode(center,i+1) )
		  {
			  number = i+1;
			  continue;
		  }
	  }
	  return number;
  }


  template< class FieldVector >
  int fixToElectrodeNumber_node (const vector<int> &corners)
  {
	  int number = 0;
	  for ( unsigned int i=0 ; i<electrode_positions_.size() ; ++i)
	  {
		  if ( fixToElectrode_node(corners,i+1) )
		  {
			  number = i+1;
			  continue;
		  }
	  }
	  return number;
  }


  // find if point is on electrode
  template< class FieldVector >
  bool onElectrode ( const FieldVector &center, int electrodenumber ) const
  {
#if 1
    // normal circular electrodes
	  const double diameter = electrode_diameter_[electrodenumber-1];
	  double electroderadius = diameter/2000;
    FieldVector electrodecenter(center);

	  electrodecenter[0]=electrode_positions_[electrodenumber-1][0];
	  electrodecenter[1]=electrode_positions_[electrodenumber-1][1];
	  electrodecenter[2]=electrode_positions_[electrodenumber-1][2];

    electrodecenter -= center;
    if( electrodecenter.two_norm() < electroderadius )
    {
      return 1;
    }
    else
      return 0;
#else
    // electrodes on unitcube
    if( electrodenumber == 1 )
	  {
      if( center[0] < 0.000000000001 )
    	  return 1;
      else
    	  return 0;
	  }
	  else if( electrodenumber == 2)
	  {
	    if( center[0] > 0.999999999999 )
		    return 1;
	    else
		    return 0;
	  }
	  else
		  return 0;
#endif
  }


  // find if part of electrode based on the corner nodes
  bool onElectrode_node ( const vector<int> &corners, int electrodenumber ) const
  {
    // corners are already sorted by default, electrode_nodes_ have to be sorted as well

    // check if all three corners are part of electrode
    int number_of_corners = 0;
	  for ( unsigned int i=0; i<electrode_positions_[electrodenumber-1].size(); i++)
    {
      if ( corners[number_of_corners] == int(electrode_positions_[electrodenumber-1][i]) )
        number_of_corners++;
      if (number_of_corners == 3)
        return 1;
    }
    
    return 0;
  }


  // find if point is on electrode
  template< class FieldVector >
  bool fixToElectrode ( const FieldVector &center, int electrodenumber )
  {
#if 1
    // normal circular electrodes
	  const double diameter = electrode_diameter_[electrodenumber-1];
	  double electroderadius = diameter/2000;
    FieldVector electrodecenter(center);

	  electrodecenter[0]=electrode_positions_[electrodenumber-1][0];
	  electrodecenter[1]=electrode_positions_[electrodenumber-1][1];
	  electrodecenter[2]=electrode_positions_[electrodenumber-1][2];

    electrodecenter -= center;
    if( electrodecenter.two_norm() < 2*electroderadius )
    {
      return 1;
    }
    else
      return 0;
#else
    // electrodes on unit cube
    if( electrodenumber == 1 )
	  {
      if( center[0] < 0.2 )
    	  return 1;
      else
    	  return 0;
	  }
	  else if( electrodenumber == 2)
	  {
	    if( center[0] > 0.8 )
		    return 1;
	    else
		    return 0;
	  }
	  else
		  return 0;
#endif
  }


  // find if part of electrode based on the corner nodes
  bool fixToElectrode_node ( const vector<int> &corners, int electrodenumber ) const
  {
    // sort corners (electrode_nodes_ are also sorted)
    vector<int> sorted_corners(corners);
    sort(sorted_corners.begin(), sorted_corners.end());

    // simply check if one of the corners is part of electrode
    int shift = 0;
	  for ( unsigned int i=0; i<electrode_positions_[electrodenumber-1].size(); ++i)
    {
      if ( sorted_corners[shift] == int(electrode_positions_[electrodenumber-1][i]) )
        return 1;

      // jump over this element corner, because it is not part of electrode
      if ( sorted_corners[shift] < int(electrode_positions_[electrodenumber-1][i]) )
      {
        if (shift==3)
          return 0;
        shift++;
        i--;
      }
    }
    
    return 0;
  }


  template< class FieldVector >
  int whichProcessShouldItBeSir ( const FieldVector &center, int numProcs)
  {
	  if (electrode_process_.size() != electrode_positions_.size())
	    writeElectrodeProcesses(numProcs, electrode_positions_.size());

	  int number = 0;
	  for ( unsigned int i=0 ; i<electrode_positions_.size() ; ++i)
	  {
	    if ( fixToElectrode(center,i+1) )
	    {
		    number = i+1;
		    continue;
	    }
	  }

	  int process = 0;
	  if (number)
	    process = electrode_process_[number-1]+1;

	  return process;
  }


  template< class FieldVector >
  int whichProcessShouldItBeSir2 ( const FieldVector &center, int numProcs, int myRank)
  {
	  if (electrode_process_.size() != electrode_positions_.size())
	    electrode_process_.resize(electrode_positions_.size());

	  int number = 0;
	  for ( unsigned int i=0 ; i<electrode_positions_.size() ; ++i)
	  {
	    if ( fixToElectrode(center,i+1) )
	    {
		    number = i+1;
		    continue;
	    }
	  }

	  if (number)
	  {
	    if (electrode_process_[number-1]==0)
	    	electrode_process_[number-1]=myRank+1;
	  }

	  return number;
  }



  // the same for electrode assignment using contributing nodes
  int whichProcessShouldItBeSir2_node ( const vector<int> &corners, int numProcs, int myRank)
  {
	  if (electrode_process_.size() != electrode_positions_.size())
	    electrode_process_.resize(electrode_positions_.size());

	  int number = 0;
	  for ( unsigned int i=0 ; i<electrode_positions_.size() ; ++i)
	  {
	    if ( fixToElectrode_node(corners,i+1) )
	    {
		    number = i+1;
		    continue;
	    }
	  }

	  if (number)
	  {
	    if (electrode_process_[number-1]==0)
	    	electrode_process_[number-1]=myRank+1;
	  }

	  return number;
  }



  void writeElectrodeProcesses( int numProcs, int numElectrodes)
  {
	  // Use half of the processes for the electrodes, so that half of the processes can work on the lower part of the mesh
	  // For the upper half of the mesh we distribute the electrodes quadratically according to x and y coordinates.
	  double min_x=electrode_positions_[0][0], min_y=electrode_positions_[0][1], min_z=electrode_positions_[0][2];
	  double max_x=electrode_positions_[0][0], max_y=electrode_positions_[0][1], max_z=electrode_positions_[0][2];
	  // Iterate through all electrodes to get minima and maxima
	  for (int i = 0; i<numElectrodes; i++)
	  {
	    if (numProcs == 1)
		  electrode_process_.push_back(0);
	    if (electrode_positions_[i][0] < min_x)
		  min_x = electrode_positions_[i][0];
	    if (electrode_positions_[i][0] > max_x)
		  max_x = electrode_positions_[i][0];
	    if (electrode_positions_[i][1] < min_y)
		  min_y = electrode_positions_[i][1];
	    if (electrode_positions_[i][1] > max_y)
		  max_y = electrode_positions_[i][1];
	    if (electrode_positions_[i][2] < min_z)
		  min_z = electrode_positions_[i][2];
	    if (electrode_positions_[i][2] > max_z)
		  max_z = electrode_positions_[i][2];
	  }

	  if (numProcs == 1)
	    return;

	  double number_of_processes_for_electrodes = (double) ceil((double) numProcs/2);
	  double sqrt_of_above = ceil(sqrt(number_of_processes_for_electrodes));
	  double interval_x = (max_x - min_x) / sqrt_of_above;
	  double interval_y = (max_y - min_y) / sqrt_of_above;

	  // Iterate through all electrodes to assign them to a process
	  for (int i = 0; i<numElectrodes; i++)
	  {
	    int p = (int) floor((electrode_positions_[i][0]-min_x) / interval_x) + sqrt_of_above*floor((electrode_positions_[i][1]-min_y) / interval_y);
	    if (p >= number_of_processes_for_electrodes)
		  p = p % (int) number_of_processes_for_electrodes;
	    electrode_process_.push_back(p);
	  }
  }


  double contactImpedance( int electrodenumber ) const
  {
    return contactimpedance_[electrodenumber-1]; // needs to be multiplied by electrode surface.
  }

  vector<int> electrodeProcessList()
  {
	  return electrode_process_;
  }


  vector<vector<double> > electrode_positions_;
  vector<int> electrode_process_;
  vector<double> electrode_diameter_;
  vector<double> contactimpedance_;
  int error_;
};

#endif // #ifndef ELECTRODE_HELPERS_HH
