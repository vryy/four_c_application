/*
LICENSE: see four_c_application/LICENSE.txt
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25/09/2024 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"

// Project includes
#include "four_c_model.h"


namespace Kratos
{

FourCModel::FourCModel(const MPI_Comm& comm, const int dim)
: mDimension(dim)
{
    mpComm = Teuchos::rcp(new Epetra_MpiComm(comm));
}

void FourCModel::CreateDiscretization(const std::string& name)
{
    auto it = mDiscretizationData.find(name);
    if (it != mDiscretizationData.end())
        KRATOS_ERROR << "4C Discretization " << name << " existed";

    mDiscretizationData[name].disc = Teuchos::rcp(new FourC::Core::FE::Discretization(name, mpComm, mDimension));
}

void FourCModel::CreateNode(const std::string& dis_name, const IndexType id, const double x, const double y, const double z)
{
    const int rank = mpComm->MyPID();
    Teuchos::RCP<FourC::Core::Nodes::Node> p_node = Teuchos::rcp(new FourC::Core::Nodes::Node(id, {x, y, z}, rank));
    auto pdisc = pGetDiscretization(dis_name);
    if (pdisc->filled())
        KRATOS_ERROR << "Discretization " << dis_name << " has been filled. No node can be added.";
    pdisc->add_node(p_node);
}

void FourCModel::CreateElement(const std::string& dis_name, const std::string& element_type, const IndexType id, const std::vector<IndexType>& nodes)
{
    if (element_type == "SOLID")
    {
        const int rank = mpComm->MyPID();
        auto p_element = Teuchos::rcp(new FourC::Discret::ELEMENTS::Solid(id, rank));
        p_element->set_node_ids(nodes.size(), nodes.data());
        auto pdisc = pGetDiscretization(dis_name);
        if (pdisc->filled())
            KRATOS_ERROR << "Discretization " << dis_name << " has been filled. No element can be added.";
        pdisc->add_element(p_element);
    }
    else
        KRATOS_ERROR << "Element type " << element_type << " is not supported";
}

void FourCModel::AddDiscretization(Teuchos::RCP<FourC::Core::FE::Discretization> pdisc)
{
    auto it = mDiscretizationData.find(pdisc->name());
    if (it != mDiscretizationData.end())
        KRATOS_ERROR << "4C Discretization " << pdisc->name() << " existed";

    mDiscretizationData[pdisc->name()].disc = pdisc;
}

void FourCModel::FillCompleteDiscretization()
{
    for (auto it = mDiscretizationData.begin(); it != mDiscretizationData.end(); ++it)
    {
        it->second.disc->fill_complete();

        it->second.mat1 = Teuchos::rcp(new FourC::Core::LinAlg::SparseMatrix(*it->second.disc->dof_row_map(), 81, true, true));
        it->second.mat2 = Teuchos::rcp(new FourC::Core::LinAlg::SparseMatrix(*it->second.disc->dof_row_map(), 81, true, true));

        it->second.vec1 = FourC::Core::LinAlg::create_vector(*it->second.disc->dof_row_map(), true);
        it->second.vec2 = FourC::Core::LinAlg::create_vector(*it->second.disc->dof_row_map(), true);
        it->second.vec3 = FourC::Core::LinAlg::create_vector(*it->second.disc->dof_row_map(), true);
    }
}

void FourCModel::EvaluateSystem(Teuchos::ParameterList& params, const std::string& dis_name)
{
    auto dt = GetDiscretizationData(dis_name);
    dt.disc->evaluate(params, dt.mat1, Teuchos::null, dt.vec1, Teuchos::null, Teuchos::null);
}

void FourCModel::EvaluateSystem(const std::string& dis_name)
{
    Teuchos::ParameterList p; // TODO fill the parameter list
    EvaluateSystem(p, dis_name);
}

} // namespace Kratos
