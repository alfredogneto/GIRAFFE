#include "Verlet.h"

#include"Database.h"

//Variáveis globais
extern
Database db;

Verlet::Verlet()
{
	neighborhood_PP = NULL;
	neighborhood_PB = NULL;
	neighborhood_BOBO = NULL;
	neighborhood_PBO = NULL;
	
	acc_steps = 0;
	n_sampling = 0;

	combine_linked_cells = false;

	//Default values for parameters of the method
	MAX_SAMPLE_RATE = 1;
	f1 = 1.0;
	f2 = 1.0;
}

Verlet::~Verlet()
{
	FreeTables();
	if (neighborhood_PP != NULL)
		delete[]neighborhood_PP;
	if (neighborhood_PB != NULL)
		delete[]neighborhood_PB;
	if (neighborhood_BOBO != NULL)
		delete[]neighborhood_BOBO;
	if (neighborhood_PBO != NULL)
		delete[]neighborhood_PBO;
}

bool Verlet::Read(FILE *f)
{
	char s[1000];

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	/*fscanf(f, "%s", s);
	if (!strcmp(s, "F1"))
	{
		fscanf(f, "%s", s);
		f1 = (float)atof(s);
	}
	else
		fsetpos(f, &pos);
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "F2"))
	{
		fscanf(f, "%s", s);
		f2 = (float)atof(s);
	}
	else
		fsetpos(f, &pos);
	fgetpos(f, &pos);*/
	fscanf(f, "%s", s);
	if (!strcmp(s, "MaxSample"))
	{
		fscanf(f, "%s", s);
		MAX_SAMPLE_RATE = atoi(s);
	}
	else
		fsetpos(f, &pos);

	return true;
}
void Verlet::Write(FILE *f)
{
	fprintf(f, "F1\t%.6e\tF2\t%.6e\tMaxSample\t%d\t", f1, f2, MAX_SAMPLE_RATE);
}

void Verlet::AllocTables()
{
	neighborhood_PP = new vector<std::array<int, 3>>[db.number_particles];
	neighborhood_PB = new vector<std::array<int, 3>>[db.number_particles];
	neighborhood_BOBO = new vector<std::array<int, 3>>[db.number_body_geometries];
	neighborhood_PBO = new vector<std::array<int, 3>>[db.number_particles];
	for (int i = 0; i < db.number_particles; i++)
	{
		neighborhood_PP[i].clear();
		neighborhood_PB[i].clear();
		neighborhood_PBO[i].clear();
	}
	for (int i = 0; i < db.number_body_geometries; i++)
	{
		neighborhood_BOBO[i].clear();
	}
}

void Verlet::FreeTables()
{
	if (neighborhood_PP != NULL)
	{
		for (int i = 0; i < db.number_particles; i++)
			neighborhood_PP[i].clear();
	}
	if (neighborhood_PB != NULL)
	{
		for (int i = 0; i < db.number_particles; i++)
			neighborhood_PB[i].clear();
	}
	if (neighborhood_BOBO != NULL)
	{
		for (int i = 0; i < db.number_body_geometries; i++)
			neighborhood_BOBO[i].clear();
	}
	if (neighborhood_PBO != NULL)
	{
		for (int i = 0; i < db.number_particles; i++)
			neighborhood_PBO[i].clear();
	}
}
void Verlet::RefreshVerletTables()
{
	FreeTables();

	//Samling rate variables - setting the value as MAX_SAMPLE_RATE
	int final_sample = MAX_SAMPLE_RATE;
	float inf = db.gcs->inc_len_factor;


	if (combine_linked_cells == false)
	{
		//Loops to update Verlet tables, according to cut-off radii and proximity between entities
//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic) 
			for (int i = 0; i < db.number_particles; i++)
			{
				for (int j = i + 1; j < db.number_particles; j++)
				{
					bool bool_near;
					int Nsampling;
					
					ComputeVerletUpdateItems(db.particles[i]->bv->x_center, db.particles[j]->bv->x_center, db.particles[i]->bv->size, db.particles[j]->bv->size,
											db.particles[i]->bv->last_center_change, db.particles[j]->bv->last_center_change, &Nsampling, &bool_near);
					//If near contact, continues
					if (bool_near)
					{
						//Searching subdivisions of particles
						for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
						{
							for (int jj = 0; jj < db.particles[j]->n_sub_bv; jj++)
							{
								ComputeVerletUpdateItems(db.particles[i]->sub_bv[ii]->x_center, db.particles[j]->sub_bv[jj]->x_center, db.particles[i]->sub_bv[ii]->size, db.particles[j]->sub_bv[jj]->size,
									db.particles[i]->sub_bv[ii]->last_center_change, db.particles[j]->sub_bv[jj]->last_center_change, &Nsampling, &bool_near);
								//If near contact, marks in Verlet table
								if (bool_near)
									neighborhood_PP[i].push_back({ j,ii,jj });
								else
								{
									if (final_sample > Nsampling)
										final_sample = Nsampling;
								}
							}
						}
					}
					else
					{
						if (final_sample > Nsampling)
							final_sample = Nsampling;
					}
				}
			}
		}

		//Loops to update Verlet tables, according to cut-off radii and proximity between entities
//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic) 
			for (int i = 0; i < db.number_particles; i++)
			{
				for (int j = 0; j < db.number_boundaries; j++)
				{
					bool bool_near;
					int Nsampling;
					ComputeVerletUpdateItems(db.particles[i]->bv->x_center, db.boundaries[j]->bv->x_center, db.particles[i]->bv->size, db.boundaries[j]->bv->size,
						db.particles[i]->bv->last_center_change, db.boundaries[j]->bv->last_center_change, &Nsampling, &bool_near);
					//If near contact, continues
					if (bool_near)
					{
						//Searching subdivisions of particles
						for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
						{
							for (int jj = 0; jj < db.boundaries[j]->n_sub_bv; jj++)
							{
								ComputeVerletUpdateItems(db.particles[i]->sub_bv[ii]->x_center, db.boundaries[j]->sub_bv[jj]->x_center, db.particles[i]->sub_bv[ii]->size, db.boundaries[j]->sub_bv[jj]->size,
									db.particles[i]->sub_bv[ii]->last_center_change, db.boundaries[j]->sub_bv[jj]->last_center_change, &Nsampling, &bool_near);
								//If near contact, marks in Verlet table
								if (bool_near)
									neighborhood_PB[i].push_back({ j,ii,jj });
								else
								{
									if (final_sample > Nsampling)
										final_sample = Nsampling;
								}
							}
						}
					}
					else
					{
						if (final_sample > Nsampling)
							final_sample = Nsampling;
					}
				}
			}
		}
		//Loops to update Verlet tables, according to cut-off radii and proximity between entities
//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic) 
			for (int i = 0; i < db.number_body_geometries; i++)
			{
				for (int j = i + 1; j < db.number_body_geometries; j++)
				{
					bool bool_near;
					int Nsampling;
					ComputeVerletUpdateItems(db.body_geometries[i]->bv->x_center, db.body_geometries[j]->bv->x_center, db.body_geometries[i]->bv->size, db.body_geometries[j]->bv->size,
						db.body_geometries[i]->bv->last_center_change, db.body_geometries[j]->bv->last_center_change, &Nsampling, &bool_near);
					//If near contact, continues
					if (bool_near)
					{
						//Searching subdivisions of body geometries
						for (int ii = 0; ii < db.body_geometries[i]->n_items; ii++)
						{
							for (int jj = 0; jj < db.body_geometries[j]->n_items; jj++)
							{
								ComputeVerletUpdateItems(db.body_geometries[i]->ptr_geom[ii]->bv->x_center, db.body_geometries[j]->ptr_geom[jj]->bv->x_center, db.body_geometries[i]->ptr_geom[ii]->bv->size, db.body_geometries[j]->ptr_geom[jj]->bv->size,
									db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change, db.body_geometries[j]->ptr_geom[jj]->bv->last_center_change, &Nsampling, &bool_near);
								//If near contact, marks in Verlet table
								if (bool_near)
									neighborhood_BOBO[i].push_back({ j,ii,jj });
								else
								{
									if (final_sample > Nsampling)
										final_sample = Nsampling;
								}
							}
						}
					}
					else
					{
						if (final_sample > Nsampling)
							final_sample = Nsampling;
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Loops to update Verlet tables, according to cut-off radii and proximity between entities
//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic) 
			for (int i = 0; i < db.number_particles; i++)
			{
				for (int j = 0; j < db.number_body_geometries; j++)
				{
					bool bool_near;
					int Nsampling;
					ComputeVerletUpdateItems(db.particles[i]->bv->x_center, db.body_geometries[j]->bv->x_center, db.particles[i]->bv->size, db.body_geometries[j]->bv->size,
						db.particles[i]->bv->last_center_change, db.body_geometries[j]->bv->last_center_change, &Nsampling, &bool_near);
					//If near contact, continues
					if (bool_near)
					{
						//Searching subdivisions of particles
						for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
						{
							for (int jj = 0; jj < db.body_geometries[j]->n_items; jj++)
							{
								ComputeVerletUpdateItems(db.particles[i]->sub_bv[ii]->x_center, db.body_geometries[j]->ptr_geom[jj]->bv->x_center, db.particles[i]->sub_bv[ii]->size, db.body_geometries[j]->ptr_geom[jj]->bv->size,
									db.particles[i]->sub_bv[ii]->last_center_change, db.body_geometries[j]->ptr_geom[jj]->bv->last_center_change, &Nsampling, &bool_near);
								//If near contact, marks in Verlet table
								if (bool_near)
									neighborhood_PBO[i].push_back({ j,ii,jj });
								else
								{
									if (final_sample > Nsampling)
										final_sample = Nsampling;
								}
							}
						}
					}
					else
					{
						if (final_sample > Nsampling)
							final_sample = Nsampling;
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}
	else
	{
//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic)
			//LinkedCells algorithm - PP contact search
			for (int i = 0; i < db.number_particles; i++)
			{
				//Pointers defined inside the loop to avoid paralellization problems
				BoundingVolume* ptr_bv;
				BoundingVolume* ptr_subbv;
				//Cell of the particle BV
				int i_loc = db.particles[i]->bv->i_loc;
				int j_loc = db.particles[i]->bv->j_loc;
				int k_loc = db.particles[i]->bv->k_loc;
				//Loops to run on 26 neighbor cells + own central cell
				for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
				{
					for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
					{
						for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
						{
							ptr_bv = db.gcs->cells_P.GetCellListRoot(cell_i, cell_j, cell_k);
							for (int ob = 0; ob < db.gcs->cells_P.GetNObjects(cell_i, cell_j, cell_k); ob++)
							{
								//Considers only if the ID of the found particle is larger than the current one
								//This avoids self-contact occurrence and doubling contact interactions
								if (ptr_bv->associated_ID > db.particles[i]->number)
								{
									int j = ptr_bv->associated_ID - 1;
									bool bool_near;
									int Nsampling;
									ComputeVerletUpdateItems(db.particles[i]->bv->x_center, db.particles[j]->bv->x_center, db.particles[i]->bv->size, db.particles[j]->bv->size,
										db.particles[i]->bv->last_center_change, db.particles[j]->bv->last_center_change, &Nsampling, &bool_near);
									//If near contact, continues
									if (bool_near)
									{
										//Searching subdivisions between particles:
										for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
										{
											//Cell of the particle subBV
											int i_subloc = db.particles[i]->sub_bv[ii]->i_loc;
											int j_subloc = db.particles[i]->sub_bv[ii]->j_loc;
											int k_subloc = db.particles[i]->sub_bv[ii]->k_loc;
											//Loops to run on 26 neighbor cells + own central cell
											for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
											{
												for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
												{
													for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
													{
														ptr_subbv = db.gcs->cells_P_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);

														for (int subob = 0; subob < db.gcs->cells_P_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
														{
															int jj = ptr_subbv->associated_sub_ID - 1;
															ComputeVerletUpdateItems(db.particles[i]->sub_bv[ii]->x_center, db.particles[j]->sub_bv[jj]->x_center, db.particles[i]->sub_bv[ii]->size, db.particles[j]->sub_bv[jj]->size,
																db.particles[i]->sub_bv[ii]->last_center_change, db.particles[j]->sub_bv[jj]->last_center_change, &Nsampling, &bool_near);
															//If near contact, marks in Verlet table
															if (bool_near)
																neighborhood_PP[i].push_back({ j,ii,jj });
															else
															{
																if (final_sample > Nsampling)
																	final_sample = Nsampling;
															}

															//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
															ptr_subbv = ptr_subbv->next;
														}
													}
												}
											}
										}
									}
									else
									{
										if (final_sample > Nsampling)
											final_sample = Nsampling;
									}
								}
								//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
								ptr_bv = ptr_bv->next;
							}
						}
					}
				}
			}
		}

//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic) 
			//LinkedCells algorithm - PB contact search
			for (int i = 0; i < db.number_particles; i++)
			{
				//Pointers defined inside the loop to avoid paralellization problems
				BoundingVolume* ptr_bv;
				BoundingVolume* ptr_subbv;
				//Cell of the particle BV
				int i_loc = db.particles[i]->bv->i_loc;
				int j_loc = db.particles[i]->bv->j_loc;
				int k_loc = db.particles[i]->bv->k_loc;
				//Loops to run on 26 neighbor cells + own central cell
				for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
				{
					for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
					{
						for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
						{
							ptr_bv = db.gcs->cells_B.GetCellListRoot(cell_i, cell_j, cell_k);
							for (int ob = 0; ob < db.gcs->cells_B.GetNObjects(cell_i, cell_j, cell_k); ob++)
							{
								int j = ptr_bv->associated_ID - 1;
								bool bool_near;
								int Nsampling;
								ComputeVerletUpdateItems(db.particles[i]->bv->x_center, db.boundaries[j]->bv->x_center, db.particles[i]->bv->size, db.boundaries[j]->bv->size,
									db.particles[i]->bv->last_center_change, db.boundaries[j]->bv->last_center_change, &Nsampling, &bool_near);
								//If near contact, continues
								if (bool_near)
								{
									//Searching subdivisions between particles/boundaries:
									for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
									{
										//Cell of the particle subBV 
										int i_subloc = db.particles[i]->sub_bv[ii]->i_loc;
										int j_subloc = db.particles[i]->sub_bv[ii]->j_loc;
										int k_subloc = db.particles[i]->sub_bv[ii]->k_loc;
										//Loops to run on 26 neighbor cells + own central cell
										for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
										{
											for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
											{
												for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
												{
													ptr_subbv = db.gcs->cells_B_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);
													for (int subob = 0; subob < db.gcs->cells_B_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
													{
														int jj = ptr_subbv->associated_sub_ID - 1;
														ComputeVerletUpdateItems(db.particles[i]->sub_bv[ii]->x_center, db.boundaries[j]->sub_bv[jj]->x_center, db.particles[i]->sub_bv[ii]->size, db.boundaries[j]->sub_bv[jj]->size,
															db.particles[i]->sub_bv[ii]->last_center_change, db.boundaries[j]->sub_bv[jj]->last_center_change, &Nsampling, &bool_near);
														//If near contact, marks in Verlet table
														if (bool_near)
															neighborhood_PB[i].push_back({ j,ii,jj });
														else
														{
															if (final_sample > Nsampling)
																final_sample = Nsampling;
														}

														//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
														ptr_subbv = ptr_subbv->next;
													}
												}
											}
										}
									}
								}
								else
								{
									if (final_sample > Nsampling)
										final_sample = Nsampling;
								}
								//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
								ptr_bv = ptr_bv->next;
							}
						}
					}
				}
			}
		}

//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic) 
			//LinkedCells algorithm - BOBO contact search
			for (int i = 0; i < db.number_body_geometries; i++)
			{
				//Pointers defined inside the loop to avoid paralellization problems
				BoundingVolume* ptr_bv;
				BoundingVolume* ptr_subbv;
				//Cell of the particle BV
				int i_loc = db.body_geometries[i]->bv->i_loc;
				int j_loc = db.body_geometries[i]->bv->j_loc;
				int k_loc = db.body_geometries[i]->bv->k_loc;
				//Loops to run on 26 neighbor cells + own central cell
				for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
				{
					for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
					{
						for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
						{
							ptr_bv = db.gcs->cells_BO.GetCellListRoot(cell_i, cell_j, cell_k);
							for (int ob = 0; ob < db.gcs->cells_BO.GetNObjects(cell_i, cell_j, cell_k); ob++)
							{
								//Considers only if the ID of the found body geometry is larger than the current one
								//This avoids self-contact occurrence and doubling contact interactions
								if (ptr_bv->associated_ID > db.body_geometries[i]->number)
								{
									int j = ptr_bv->associated_ID - 1;
									bool bool_near;
									int Nsampling;
									ComputeVerletUpdateItems(db.body_geometries[i]->bv->x_center, db.body_geometries[j]->bv->x_center, db.body_geometries[i]->bv->size, db.body_geometries[j]->bv->size,
										db.body_geometries[i]->bv->last_center_change, db.body_geometries[j]->bv->last_center_change, &Nsampling, &bool_near);
									//If near contact, continues
									if (bool_near)
									{
										//Searching subdivisions between body geometries:
										for (int ii = 0; ii < db.body_geometries[i]->n_items; ii++)
										{
											//Cell of the particle subBV
											int i_subloc = db.body_geometries[i]->ptr_geom[ii]->bv->i_loc;
											int j_subloc = db.body_geometries[i]->ptr_geom[ii]->bv->j_loc;
											int k_subloc = db.body_geometries[i]->ptr_geom[ii]->bv->k_loc;
											//Loops to run on 26 neighbor cells + own central cell
											for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
											{
												for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
												{
													for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
													{
														ptr_subbv = db.gcs->cells_BO_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);

														for (int subob = 0; subob < db.gcs->cells_BO_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
														{
															int jj = ptr_subbv->associated_sub_ID - 1;
															ComputeVerletUpdateItems(db.body_geometries[i]->ptr_geom[ii]->bv->x_center, db.body_geometries[j]->ptr_geom[jj]->bv->x_center, db.body_geometries[i]->ptr_geom[ii]->bv->size, db.body_geometries[j]->ptr_geom[jj]->bv->size,
																db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change, db.body_geometries[j]->ptr_geom[jj]->bv->last_center_change, &Nsampling, &bool_near);
															//If near contact, marks in Verlet table
															if (bool_near)
																neighborhood_BOBO[i].push_back({ j,ii,jj });
															else
															{
																if (final_sample > Nsampling)
																	final_sample = Nsampling;
															}

															//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
															ptr_subbv = ptr_subbv->next;
														}
													}
												}
											}
										}
									}
									else
									{
										if (final_sample > Nsampling)
											final_sample = Nsampling;
									}
								}

								//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
								ptr_bv = ptr_bv->next;
							}
						}
					}
				}
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#pragma omp parallel
		{
//#pragma omp for schedule(dynamic) 
			//LinkedCells algorithm - PBO contact search
			for (int i = 0; i < db.number_particles; i++)
			{
				//Pointers defined inside the loop to avoid paralellization problems
				BoundingVolume* ptr_bv;
				BoundingVolume* ptr_subbv;
				//Cell of the particle BV
				int i_loc = db.particles[i]->bv->i_loc;
				int j_loc = db.particles[i]->bv->j_loc;
				int k_loc = db.particles[i]->bv->k_loc;
				//Loops to run on 26 neighbor cells + own central cell
				for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
				{
					for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
					{
						for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
						{
							ptr_bv = db.gcs->cells_BO.GetCellListRoot(cell_i, cell_j, cell_k);
							for (int ob = 0; ob < db.gcs->cells_BO.GetNObjects(cell_i, cell_j, cell_k); ob++)
							{
								int j = ptr_bv->associated_ID - 1;
								bool bool_near;
								int Nsampling;
								ComputeVerletUpdateItems(db.particles[i]->bv->x_center, db.body_geometries[j]->bv->x_center, db.particles[i]->bv->size, db.body_geometries[j]->bv->size,
									db.particles[i]->bv->last_center_change, db.body_geometries[j]->bv->last_center_change, &Nsampling, &bool_near);
								//If near contact, continues
								if (bool_near)
								{
									//Searching subdivisions between particles/body_geometries:
									for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
									{
										//Cell of the particle subBV 
										int i_subloc = db.particles[i]->sub_bv[ii]->i_loc;
										int j_subloc = db.particles[i]->sub_bv[ii]->j_loc;
										int k_subloc = db.particles[i]->sub_bv[ii]->k_loc;
										//Loops to run on 26 neighbor cells + own central cell
										for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
										{
											for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
											{
												for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
												{
													ptr_subbv = db.gcs->cells_BO_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);
													for (int subob = 0; subob < db.gcs->cells_BO_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
													{
														int jj = ptr_subbv->associated_sub_ID - 1;
														ComputeVerletUpdateItems(db.particles[i]->sub_bv[ii]->x_center, db.body_geometries[j]->ptr_geom[jj]->bv->x_center, db.particles[i]->sub_bv[ii]->size, db.body_geometries[j]->ptr_geom[jj]->bv->size,
															db.particles[i]->sub_bv[ii]->last_center_change, db.body_geometries[j]->ptr_geom[jj]->bv->last_center_change, &Nsampling, &bool_near);
														//If near contact, marks in Verlet table
														if (bool_near)
															neighborhood_PBO[i].push_back({ j,ii,jj });
														else
														{
															if (final_sample > Nsampling)
																final_sample = Nsampling;
														}

														//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
														ptr_subbv = ptr_subbv->next;
													}
												}
											}
										}
									}
								}
								else
								{
									if (final_sample > Nsampling)
										final_sample = Nsampling;
								}
								//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
								ptr_bv = ptr_bv->next;
							}
						}
					}
				}
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////
	}
	n_sampling = final_sample;
}

void Verlet::SetSamplingRate()
{
	//int final_sample = INT_MAX;
	//float inf = db.gcs->inc_len_factor;
	//
	//if (MAX_SAMPLE_RATE != 1)
	//{
	//	//Searching particles
	//	for (int i = 0; i < db.number_particles; i++)
	//	{
	//		float last_change1 = sqrt(db.particles[i]->bv->last_center_change[0] * db.particles[i]->bv->last_center_change[0] +
	//			db.particles[i]->bv->last_center_change[1] * db.particles[i]->bv->last_center_change[1] +
	//			db.particles[i]->bv->last_center_change[2] * db.particles[i]->bv->last_center_change[2]);
	//		if (last_change1 != 0.0f)
	//		{
	//			int temp_sample1 = (int)((inf/(2.0f*(1.0f+inf)))*db.particles[i]->bv->size / last_change1);
	//			if (temp_sample1 < final_sample)
	//				final_sample = temp_sample1;
	//		}
	//		//Searching sub-particles
	//		for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
	//		{
	//			float last_change2 = sqrt(db.particles[i]->sub_bv[ii]->last_center_change[0] * db.particles[i]->sub_bv[ii]->last_center_change[0] +
	//				db.particles[i]->sub_bv[ii]->last_center_change[1] * db.particles[i]->sub_bv[ii]->last_center_change[1] +
	//				db.particles[i]->sub_bv[ii]->last_center_change[2] * db.particles[i]->sub_bv[ii]->last_center_change[2]);
	//			if (last_change2 != 0.0f)
	//			{
	//				int temp_sample2 = (int)((inf / (2.0f*(1.0f + inf)))*db.particles[i]->sub_bv[ii]->size / last_change2);
	//				if (temp_sample2 < final_sample)
	//					final_sample = temp_sample2;
	//			}
	//		}
	//	}
	//	//Searching boundaries
	//	for (int i = 0; i < db.number_boundaries; i++)
	//	{
	//		float last_change1 = sqrt(db.boundaries[i]->bv->last_center_change[0] * db.boundaries[i]->bv->last_center_change[0] +
	//			db.boundaries[i]->bv->last_center_change[1] * db.boundaries[i]->bv->last_center_change[1] +
	//			db.boundaries[i]->bv->last_center_change[2] * db.boundaries[i]->bv->last_center_change[2]);
	//		if (last_change1 != 0.0f)
	//		{
	//			int temp_sample1 = (int)((inf / (2.0f*(1.0f + inf)))*db.boundaries[i]->bv->size / last_change1);
	//			if (temp_sample1 < final_sample)
	//				final_sample = temp_sample1;
	//		}
	//		//Searching sub-boundaries
	//		for (int ii = 0; ii < db.boundaries[i]->n_sub_bv; ii++)
	//		{
	//			float last_change2 = sqrt(db.boundaries[i]->sub_bv[ii]->last_center_change[0] * db.boundaries[i]->sub_bv[ii]->last_center_change[0] +
	//				db.boundaries[i]->sub_bv[ii]->last_center_change[1] * db.boundaries[i]->sub_bv[ii]->last_center_change[1] +
	//				db.boundaries[i]->sub_bv[ii]->last_center_change[2] * db.boundaries[i]->sub_bv[ii]->last_center_change[2]);
	//			if (last_change2 != 0.0f)
	//			{
	//				int temp_sample2 = (int)((inf / (2.0f*(1.0f + inf)))*db.boundaries[i]->sub_bv[ii]->size / last_change2);
	//				if (temp_sample2 < final_sample)
	//					final_sample = temp_sample2;
	//			}
	//		}
	//	}
	//	//Searching body geometries
	//	for (int i = 0; i < db.number_body_geometries; i++)
	//	{
	//		float last_change1 = sqrt(db.body_geometries[i]->bv->last_center_change[0] * db.body_geometries[i]->bv->last_center_change[0] +
	//			db.body_geometries[i]->bv->last_center_change[1] * db.body_geometries[i]->bv->last_center_change[1] +
	//			db.body_geometries[i]->bv->last_center_change[2] * db.body_geometries[i]->bv->last_center_change[2]);
	//		if (last_change1 != 0.0f)
	//		{
	//			int temp_sample1 = (int)((inf / (2.0f*(1.0f + inf)))*db.body_geometries[i]->bv->size / last_change1);
	//			if (temp_sample1 < final_sample)
	//				final_sample = temp_sample1;
	//		}
	//		//Searching sub-body geometries
	//		for (int ii = 0; ii < db.body_geometries[i]->n_items; ii++)
	//		{
	//			float last_change2 = sqrt(db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change[0] * db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change[0] +
	//				db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change[1] * db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change[1] +
	//				db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change[2] * db.body_geometries[i]->ptr_geom[ii]->bv->last_center_change[2]);
	//			if (last_change2 != 0.0f)
	//			{
	//				int temp_sample2 = (int)((inf / (2.0f*(1.0f + inf)))*db.body_geometries[i]->ptr_geom[ii]->bv->size / last_change2);
	//				if (temp_sample2 < final_sample)
	//					final_sample = temp_sample2;
	//			}
	//		}
	//	}

	//	//printf("SAMPLE %d\n--------------------------------", final_sample);
	//	if (final_sample > MAX_SAMPLE_RATE)
	//		final_sample = MAX_SAMPLE_RATE;
	//	if (final_sample == 0)
	//		final_sample = 1;
	//}
	//else
	//	final_sample = 1;
	//
	//n_sampling = final_sample;
}

void Verlet::ComputeVerletUpdateItems(float* C_i, float* C_j, float size_i, float size_j, float* l_i, float* l_j, int* Nsampling, bool* bool_near)
{
	//vector pointing from center j to center i
	float v_ij[3];	
	v_ij[0] = C_i[0] - C_j[0];
	v_ij[1] = C_i[1] - C_j[1];
	v_ij[2] = C_i[2] - C_j[2];
	//center-center distance (scalar)
	float d_ij = sqrt(v_ij[0] * v_ij[0] + v_ij[1] * v_ij[1] + v_ij[2] * v_ij[2]);
	//inflated sum of BV's sizes
	float cutoff_ij = (size_i + size_j);
	float num = (d_ij - cutoff_ij);

	//if num <= 0, there is proximity enough to include the pair in the list 
	if (num <= 0.0f)
	{
		bool inside_domain;
		if (C_i[0] >= db.gcs->min_x && C_i[0] <= db.gcs->max_x &&
			C_i[1] >= db.gcs->min_y && C_i[1] <= db.gcs->max_y &&
			C_i[2] >= db.gcs->min_z && C_i[2] <= db.gcs->max_z &&
			C_j[0] >= db.gcs->min_x && C_j[0] <= db.gcs->max_x &&
			C_j[1] >= db.gcs->min_y && C_j[1] <= db.gcs->max_y &&
			C_j[2] >= db.gcs->min_z && C_j[2] <= db.gcs->max_z)
			inside_domain = true;
		else
			inside_domain = false;
		if (inside_domain)
			*bool_near = true;
		else
			*bool_near = false;
		//*Nsampling is not modified, as long as it is not going to be used in RefreshVerlet function
	}
	else
	{
		*bool_near = false;
		float den = -1.0f * ((l_i[0] - l_j[0]) * v_ij[0] + (l_i[1] - l_j[1]) * v_ij[1] + (l_i[2] - l_j[2]) * v_ij[2]) / d_ij;
		if (den > 0.0f)
		{
			//temp is another variation of 'num', employed with a distinct increase len factor. This is to avoid establishing too low Nsampling
			//if the same inc len factor is used for both inclusion on Verlet list and sampling, it is quite often to have Nsampling = 1
			//0.5f * db.gcs->inc_len_factor is a possibility. Another values could be tried, also taking some info on interface laws, etc.
			float temp = ((d_ij - cutoff_ij * (1.0f + 0.5f * db.gcs->inc_len_factor) / (1.0f + db.gcs->inc_len_factor)) / den);
			if (temp < (float)INT_MAX)
				*Nsampling = (int)temp;
			else
				*Nsampling = MAX_SAMPLE_RATE;

			if (*Nsampling == 0)
			{
				//Printing info
				/*db.myprintf("\nv_ij (%f,%f,%f)\n", v_ij[0], v_ij[1], v_ij[2]);
				db.myprintf("l_i (%f,%f,%f)\n", l_i[0], l_i[1], l_i[2]);
				db.myprintf("l_j (%f,%f,%f)\n", l_j[0], l_j[1], l_j[2]);
				db.myprintf("d_ij %f\n", d_ij);
				db.myprintf("cutoff_ij %f\n", cutoff_ij);
				db.myprintf("num %f\n", num);
				db.myprintf("den %f\n", den);
				db.myprintf("Nsampling %d\n", *Nsampling);*/
				*Nsampling = 1;

			}
			
		}
		else
			*Nsampling = MAX_SAMPLE_RATE;
		
	}
}
