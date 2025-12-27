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
: mDimension(dim), mComm(comm)
{}

FourCModel::~FourCModel()
{}

void FourCModel::CreateDiscretization(const std::string& name)
{
    auto it = mDiscretizationData.find(name);
    if (it != mDiscretizationData.end())
        KRATOS_ERROR << "4C Discretization " << name << " existed";

    mDiscretizationData[name].disc = std::make_shared<FourC::Core::FE::Discretization>(name, mComm, mDimension);
}

void FourCModel::CreateNode(const std::string& dis_name, const IndexType id, const double x, const double y, const double z)
{
    const int rank = FourC::Core::Communication::my_mpi_rank(mComm);
    std::vector<double> coords = {x, y, z};
    const auto p_node = std::make_shared<FourC::Core::Nodes::Node>(id, coords, rank);
    auto pdisc = pGetDiscretization(dis_name);
    if (pdisc->filled())
        KRATOS_ERROR << "Discretization " << dis_name << " has been filled. No node can be added.";
    std::span<const double, 3> xv(coords.data(), 3);
    pdisc->add_node(xv, id, p_node);
}

void FourCModel::CreateElement(const std::string& dis_name, const std::string& element_type, const IndexType id, const std::vector<IndexType>& nodes)
{
    const int rank = FourC::Core::Communication::my_mpi_rank(mComm);
    auto pdisc = pGetDiscretization(dis_name);
    if (pdisc->filled())
        KRATOS_ERROR << "Discretization " << dis_name << " has been filled. No element can be added.";

    std::shared_ptr<FourC::Core::Elements::Element> p_element;
    if (element_type == "SOLID")
    {
        p_element = std::make_shared<FourC::Discret::Elements::Solid>(id, rank);
    }
    else
        KRATOS_ERROR << "Element type " << element_type << " is not supported";

    p_element->set_node_ids(nodes.size(), nodes.data());
    pdisc->add_element(p_element);
}

void FourCModel::AddDiscretization(std::shared_ptr<FourC::Core::FE::Discretization> pdisc)
{
    auto it = mDiscretizationData.find(pdisc->name());
    if (it != mDiscretizationData.end())
        KRATOS_ERROR << "4C Discretization " << pdisc->name() << " existed";

    mDiscretizationData[pdisc->name()].disc = pdisc;
}

void FourCModel::FillComplete()
{
    for (auto it = mDiscretizationData.begin(); it != mDiscretizationData.end(); ++it)
    {
        it->second.disc->fill_complete();

        it->second.mat1 = std::make_shared<FourC::Core::LinAlg::SparseMatrix>(*it->second.disc->dof_row_map(), 81, true, true);
        it->second.mat2 = std::make_shared<FourC::Core::LinAlg::SparseMatrix>(*it->second.disc->dof_row_map(), 81, true, true);

        it->second.vec1 = std::make_shared<FourC::Core::LinAlg::Vector<double>>(*it->second.disc->dof_row_map(), true);
        it->second.vec2 = std::make_shared<FourC::Core::LinAlg::Vector<double>>(*it->second.disc->dof_row_map(), true);
        it->second.vec3 = std::make_shared<FourC::Core::LinAlg::Vector<double>>(*it->second.disc->dof_row_map(), true);
    }
}

void FourCModel::SetZeroState(const std::string& dis_name, const unsigned nds,
        const std::string& state_name)
{
    auto pdisc = GetDiscretizationData(dis_name).disc;
    auto zero_state = std::make_shared<FourC::Core::LinAlg::Vector<double>>(*pdisc->dof_row_map(), true);
    pdisc->set_state(nds, state_name, *zero_state);
}

void FourCModel::SetState(const std::string& dis_name, const unsigned nds,
        const std::string& state_name, std::shared_ptr<const FourC::Core::LinAlg::Vector<double> > state)
{
    auto pdisc = GetDiscretizationData(dis_name).disc;
    pdisc->set_state(nds, state_name, *state);
}

void FourCModel::Evaluate(Teuchos::ParameterList& params, const std::string& dis_name)
{
    auto dt = GetDiscretizationData(dis_name);
    dt.disc->evaluate(params, dt.mat1, dt.mat2, dt.vec1, dt.vec2, dt.vec3);
}

} // namespace Kratos
